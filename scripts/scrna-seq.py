
import os
from random import seed
import scanpy as sc
import scrublet as scr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import median_abs_deviation
import numpy as np

seed(24032026)

# User-defined parameters
PROJECT_DIR = "/work/TALC/mdsc519_2026w/students/jamie/Dragon"
CK_DIR = os.path.join(PROJECT_DIR, "data", "output", "cellranger", "SRR24952454")
POST_DIR = os.path.join(PROJECT_DIR, "data", "output", "cellranger", "SRR24952453")
CELLRANGER_OUTDIR = "/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/output/cellranger"
FIG_OUTPUTS_DIR = os.path.join(PROJECT_DIR, "output", "figures")

# Configure scanpy settings
sc.settings.figdir = os.path.join(FIG_OUTPUTS_DIR)
sc.set_figure_params(figsize=(8, 8), dpi=100)
plt.rcParams["figure.figsize"] = (8, 8)
plt.rcParams["figure.dpi"] = 100

samples = {
    "post": "SRR24952453", "ck": "SRR24952454"
}

exp_cell_lookup = {
    750: 0.004,
    1500: 0.008,
    2500: 0.016,
    3500: 0.023,
    4500: 0.031
}

# Functions for QC and filtering
def get_expected_doublet_rate(n_cells):
    # Find the closest key in the lookup dictionary
    expected_rate = 0.039
    
    if n_cells <= 750:
        expected_rate = 0.004
    elif n_cells <= 1500:
        expected_rate = 0.008
    elif n_cells <= 2500:
        expected_rate = 0.016
    elif n_cells <= 3500:
        expected_rate = 0.023
    elif n_cells <= 4500:
        expected_rate = 0.031
        
    return expected_rate

# Functions for QC and filtering
def per_sample_qc(adata, sample_name):
    print(f"\nQC metrics for sample {sample_name}:")
    print(f"  Number of cells: {adata.n_obs}")
    print(f"  Number of genes: {adata.n_vars}")
    adata.var["mt"] = adata.var_names.str.startswith(("MT-","mt-"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        log1p=True,
        inplace=True,
        percent_top=[20, 50, 100]
    )
    
    with plt.rc_context({"figure.figsize": (15, 6)}):
        sc.pl.violin(
            adata, 
            ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
            # log=True,
            jitter=0.4,
            multi_panel=True,
            stripplot=True,
            size=3,
            palette="Set2",
            save=f"_qc_{sample_name}.png"
        )
        plt.title("Total UMI Counts per Cell", fontsize=14)
        plt.ylabel("log10(Total UMI counts)", fontsize=12)
        plt.xlabel("", fontsize=12)
        
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save=f"_total_counts_vs_n_genes_{sample_name}.png")
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save=f"_total_counts_vs_pct_counts_mt_{sample_name}.png")
    
    print(f"Figures saved for sample {sample_name} QC metrics.")
    
def iqr_bounds(series, multiplier=1.5):
    # Calculate the interquartile range (IQR) and determine the lower and upper bounds for outliers based on the IQR method. The multiplier parameter determines how far from the IQR the bounds are set 
    q1 = np.percentile(series, 25)
    q3 = np.percentile(series, 75)
    iqr = q3 - q1
    lower_bound = q1 - multiplier * iqr
    upper_bound = q3 + multiplier * iqr
    return lower_bound, upper_bound


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def per_sample_filtering(adata, sample_name, method="iqr", multiplier=1.5):
    # Paper used quartile-threshold criteria, but exact cutoffs are not specified. 
    # 
    filtered= adata.copy()
    # IQR
    if method == "iqr":
        for metric in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
            lower_bound, upper_bound = iqr_bounds(filtered.obs[metric], multiplier)
            print(f"  {metric}: [{lower_bound:.3f}, {upper_bound:.3f}]")
            filtered = filtered[
                (filtered.obs[metric] >= lower_bound) &
                (filtered.obs[metric] <= upper_bound)
            ].copy()
        print(f"  Filtered number of cells: {filtered.n_obs}")
        return filtered
    # 
    # Here, we will use a more systematic approach based on MAD (Median Absolute Deviation) to identify outliers in the QC metrics and filter out low-quality cells.
    for metric in ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_in_top_20_genes"]:
        M = filtered.obs[metric]
        median = np.median(M)
        mad = median_abs_deviation(M)
        lower = median - 5 * mad
        upper = median + 5 * mad
        print(f"  {metric}: [{lower:.3f}, {upper:.3f}]")

        filtered.obs["outlier"] = (
            is_outlier(filtered, "log1p_total_counts", 5)
            | is_outlier(filtered, "log1p_n_genes_by_counts", 5)
            | is_outlier(filtered, "pct_counts_in_top_20_genes", 5)
        )
        
        print(f"  Outliers based on total counts, n_genes_by_counts, and pct_counts_in_top_20_genes: {filtered.obs.outlier.value_counts()}")
        
        filtered.obs["mt_outlier"] = is_outlier(filtered, "pct_counts_mt", 3) | (
            filtered.obs["pct_counts_mt"] > 8
        )
        print(f"  Outliers based on pct_counts_mt: {filtered.obs.mt_outlier.value_counts()}")
        
        filtered = filtered[~filtered.obs.outlier & ~filtered.obs.mt_outlier].copy()
        print(f"  Filtered number of cells: {filtered.n_obs}")
        return filtered
    
# scrublet to identify doublets - not sure if we have to do this since we already have the cellranger output, which should have already filtered out doublets. But could be worth trying to see if we can identify any additional doublets that were missed by cellranger.
def identify_doublets(adata, expected_doublet_rate=0.06, sample_name=None, umap=False):
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate, threshold=None)
    if sample_name is not None:
        sc.pl.scrublet_score_distribution(scrub, save=f"_{sample_name}.png")
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets
    print(f"Predicted doublets: {predicted_doublets.sum()} out of {adata.n_obs} cells")
    # TODO: finish the inspection and choice of doublet threshold
#     if umap:
#         tmp = adata.copy()
#         sc.pp.normalize_total(tmp, target_sum=1e4, inplace=True)
#         sc.pp.log1p(tmp)
#         sc.pp.highly_variable_genes(tmp, n_top_genes=2000)
#         tmp = tmp[:, tmp.var["highly_variable"]].copy()
#         sc.pp.scale(tmp, max_value=10)
#         sc.tl.pca(tmp)
#         sc.pp.neighbors(tmp)
#         sc.tl.umap(tmp)
        
#         tmp.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(str) # convert to string for plotting
#         sc.pl.umap(tmp, color=["predicted_doublet", "doublet_score"], save=f"doublet_viz_{sample_name}.png")

#     adata.obs["possible_doublets"] = adata.obs["predicted_doublet"].copy()
    
#     # stricter threshold for doublets
#     inspect = adata.obs["possible_doublets"].values
#     inspect_scores = adata.obs.loc[inspect, "doublet_score"]
    
#     strict_cut = np.quantile(inspect_scores, 0.5)
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    return adata
    
    
    
    
# Main analysis workflow
#########################################################################

# Redirect all output to a log file
# log_file = open(f"{OUTPUT_DIR}/output.log", "w")
# sys.stdout = log_file

adatas = []
for sample, srr_id in samples.items():
    h5_path = os.path.join(CELLRANGER_OUTDIR, sample, "outs", "filtered_feature_bc_matrix.h5")
    adata = sc.read_10x_mtx(h5_path, var_names="gene_symbols")
    adata.var_names_make_unique() # make gene names unique by adding a suffix to duplicate names, which is important for downstream analyses that require unique gene identifiers. This step ensures that each gene can be uniquely identified and prevents issues that may arise from having duplicate gene names in the dataset.
    adata.obs["sample_id"] = srr_id
    adata.obs["condition"] = sample
    adatas.append(adata)

    # for each adata, comput qc, remove outliers. Then concatenate the adatas together for downstream analyses. 
    adata = per_sample_qc(adata, sample)
    adata.raw = adata.copy() # store the original data in the raw attribute of the adata object, so that we can use it later for plotting and other analyses. This is important because the filtering step changes the values of the genes, so we want to keep a copy of the original data for reference.
    # Use more relaxed filtering criteria (3 MADs) to retain more cells for downstream analyses, as the paper does not specify the exact filtering criteria used. 
    # This allows us to keep more cells in the dataset while still removing clear outliers that may represent low-quality cells or technical artifacts. We can always adjust the multiplier parameter later if we find that we are retaining too many low-quality cells or if we want to be more stringent in our filtering.
    adata = per_sample_filtering(adata, sample, multiplier=3.0) 
    adata.obs["filtered"] = True
    
    expected_doublet_rate = get_expected_doublet_rate(adata.n_obs)
    adata = identify_doublets(adata, expected_doublet_rate=expected_doublet_rate)
    
    adatas.append(adata)

adata = sc.concat(adatas, join="outer", label="sample_id", keys=list(samples.keys()))

#
# Normalize the data
# ---------------------------------------------------------------------
print("\n\Normalizing the data")
adata.raw = adata.copy() # store the original data in the raw attribute of the adata object, so that we can use it later for plotting and other analyses. This is important because the normalization step changes the values of the genes, so we want to keep a copy of the original data for reference.
sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
sc.pp.log1p(adata)

print(str(adata.var.columns.tolist()))


# Identify highly variable genes
# ---------------------------------------------------------------------
print("\n\Identifying highly variable genes")
#https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html
sc.pp.highly_variable_genes(
    adata, # use the filtered adata?  to identify highly variable genes, as the normalization step does not change the variance of the genes, so we can use the original data to identify highly variable genes. (Get a warning when using normalized data) - flavor seurat expects normalized data, seurat_v3 expects raw counts data
    flavor="seurat",
    # span=0.5, 
    # n_top_genes=200,
    # inplace=False
    # min_mean=0.05, 
    # max_mean=1, 
    # min_disp=0.8,
    )
sc.pl.highly_variable_genes(adata, save=".png")
print(f"Number of highly variable genes: {adata.var['highly_variable'].sum()}")

# adata = adata[:, adata.var["highly_variable"]].copy() # subset the normalized adata to only include the highly variable genes, as these are the genes that will be used for downstream analyses such as dimensionality reduction and clustering. This step is important because it reduces the dimensionality? of the data and focuses on the most informative genes. A lot of cells have zero counts for many genes, so including all genes would add a lot of noise to the data and make it harder to identify meaningful patterns. By selecting only the highly variable genes, we can improve the signal-to-noise ratio and enhance the ability to detect biologically relevant clusters and patterns in the data.

# print(adata.var["highly_variable"].value_counts())
# print(adata.var[adata.var["highly_variable"]].head())
# print(f"Number of highly variable genes: {hvg.sum()}")

# Calculate means and variances of genes for plotting
# adata.var["means"] = np.array(adata.X.mean(axis=0)).flatten()
# if hasattr(adata.X, 'toarray'):  # sparse matrix
#     adata.var["variances"] = np.array(adata.X.toarray().var(axis=0)).flatten()
# else:  # dense matrix
#     adata.var["variances"] = np.array(adata.X.var(axis=0)).flatten()

# print(adata.var["means"][:5]) # print the mean expression of the first 5 genes
# print(adata.var["variances"][:5]) # print the variance of the first 5 genes
print(adata.var) # print the variable names

# plot mean vs variance of genes, highlighting the highly variable genes - DOESN'T WORK for flavor=seurat, would have to manually calculate or use seurat_v3.
# sc.pp.highly_variable_genes(
#     adata,
#     flavor="seurat_v3",
#     n_top_genes=2000
# )
# sns.scatterplot(
#     x=adata.var["means"],
#     y=adata.var["variances"],
#     hue=adata.var["highly_variable"],
#     palette={True: "red", False: "blue"},
#     alpha=0.5
# )
# plt.xlabel("Mean expression", fontsize=14)
# plt.ylabel("Variance", fontsize=14)
# plt.legend(title="Highly Variable")
# plt.title()
# plt.savefig(f"{OUTPUT_DIR}/5-mean_vs_variance.png", bbox_inches='tight')


adata = adata[:, adata.var["highly_variable"]].copy()
print(f"Shape of highly variable gene matrix: {adata.shape}")

# SCALE ?
# Remove background noise? better clustering?
# sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

sc.pp.scale(adata, max_value=10) # scale the data to unit variance and zero mean, with a maximum value of 10 to prevent extreme values from dominating the analysis. This step is important because it ensures that all genes are on the same scale and that highly expressed genes do not dominate the downstream analyses such as PCA and clustering. Scaling makes analyses more sensitive to subtle differences in gene expression across cells.


################################################################################
# Dimensionality Reduction and Clustering
################################################################################
# Run PCA and visualize the variance explained by each principal component
# ---------------------------------------------------------------------
print("\n\nRun PCA and visualize the variance explained by each principal component")
sc.pp.pca(adata, n_comps=50, svd_solver="arpack") 
# sc.tl.pca(adata, svd_solver="arpack") This version is deprecated? use sc.pp.pca instead
sc.pl.pca_variance_ratio(adata, log=True, save="6.png")

sc.pl.pca_overview(adata, save="6.png")

# Elbow plot to show variance explained by each principal component (first 40 PCs)
variance = adata.uns['pca']['variance'][:40]
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(variance)+1), variance, 'bo-', linewidth=2, markersize=6)
plt.xlabel('Principal Component', fontsize=12)
plt.ylabel('Variance Explained', fontsize=12)
plt.title('Elbow Plot', fontsize=14)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIG_OUTPUTS_DIR}/pca_elbow_plot.png", bbox_inches='tight')
plt.close()

n_PCs = 6 # choose the number of PCs to use for downstream analyses based on the elbow plot and the variance explained by each PC. In this case, we can see that the first 4 PCs explain a significant amount of variance in the data, so we will use these for clustering and visualization.

# Construct kNN graph on pca and perform clustering using the Leiden algorithm
# --------------------------------------------------------------------------------
print("\n\nConstruct kNN graph and perform clustering using the Leiden algorithm")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_PCs)
resolutions = [0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
res_table = []
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", n_iterations=2, flavor="igraph", directed=False)
    print("Available obs columns:", adata.obs.columns.tolist())
    print(f"Leiden clustering with resolution {res}:")
    print(adata.obs[f"leiden_{res}"].value_counts())
    res_table.append(
        {
            "resolution": res,
            "n_clusters": adata.obs[f"leiden_{res}"].nunique(),
            "size_of_largest_cluster": adata.obs[f"leiden_{res}"].value_counts().max(),
            "average_cluster_size": adata.obs[f"leiden_{res}"].value_counts().mean(),
            "smallest_cluster_size": adata.obs[f"leiden_{res}"].value_counts().min()
        }
    )
    
print(f"\nResolution table: \n{pd.DataFrame(res_table)}")
# Choose resolution
# -----------------------------------------------------------------------------
# choose a resolution based on the number of clusters and their sizes. A resolution of 1.0 seems to give a reasonable number of clusters (10) with a good balance between cluster sizes (largest cluster has 100 cells, average cluster size is 50 cells, and smallest cluster has 20 cells).
print("\n\nQuestion 8: Choose resolution")
chosen_res = 0.50
print(f"\nChosen resolution: {chosen_res}")
# print the number of clusters and their sizes for the chosen resolution
print(adata.obs[f"leiden_{chosen_res}"].value_counts())

# 9. UMAP embedding and plot cells
# ---------------------------------------------------------------------
#https://scanpy.readthedocs.io/en/stable/tutorials/plotting/core.html
print("\n\nQuestion 9: UMAP embedding and plot")
sc.tl.umap(adata)
# colour by cluster
# sc.pl.umap(adata, 
#         color=f"leiden_1.0", 
#     #    add_outline=True,
#         legend_loc="on data",
#         legend_fontsize=12,
#         alpha=0.5,
#     #    legend_fontoutline=2,
#         frameon=False,
#         title=f"UMAP colored by Leiden clusters (resolution=1.0)",
#         save="_clustering9.png")
# colour by total UMI counts
sc.pl.umap(adata, 
        color="total_counts", 
        frameon=False,
        save="_total_counts.png")
# colour by mitochondrial gene percentage    
sc.pl.umap(adata, color="pct_counts_mt", frameon=False, save="_pct_counts_mt.png")

# UMAP showing clusters for chosen resolution
# ---------------------------------------------------------------------
print("\n\nUMAP showing clusters for chosen resolution")
sc.pl.umap(adata, color=f"leiden_{chosen_res}", alpha=0.5, frameon=False, legend_loc="on data", legend_fontsize=12, title=f"UMAP colored by Leiden clusters (resolution={chosen_res})", save=".png")


##############################################################################
# PART 4 -  Cluster Annotation
################################################################################
# followed: https://scanpy.readthedocs.io/en/stable/tutorials/plotting/core.html
# 11. Identify marker genes for each cluster
# ---------------------------------------------------------------------
print("\n\nIdentify marker genes for each cluster")
sc.tl.dendrogram(adata, groupby=f"leiden_{chosen_res}")

# gene_marker_dict = {
#     "CD4+ T cells": ["CD3D", "CD3E", "IL7R", "CD4"],
#     "CD8+ T cells": ["CD3D", "CD3E", "CD8A", "CD8B"],
#     "NK cells": ["NKG7", "GNLY", "KLRD1"],
#     "B cells": ["MS4A1", "CD79A", "CD79B"],
#     "CD14+ Monocytes": ["CD14", "LYZ", "CST3"],
#     "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7"],
#     "Dendritic cells": ["FCER1A", "CST3"],
#     "Megakaryocytes": ["PPBP", "PF4"]
# }

# Create cluster-to-cell-type mapping based on marker gene expression

# # Add cell type annotation to adata
# adata.obs["cell_type"] = adata.obs[f"leiden_{chosen_res}"].astype(str).map(cluster_labels)

# Create marker genes list for dotplot (ordered by cluster)
# marker_genes_list = [
#     ["CD3D", "CD3E", "IL7R", "CD4"],        # CD4+ T cells
#     ["CD3D", "CD3E", "CD8A", "CD8B"],       # CD8+ T cells
#     ["NKG7", "GNLY", "KLRD1"],              # NK cells
#     ["MS4A1", "CD79A", "CD79B"],            # B cells
#     ["CD14", "LYZ", "CST3"],                # CD14+ Monocytes
#     ["FCGR3A", "MS4A7"],                    # FCGR3A+ Monocytes
#     ["FCER1A", "CST3"],                     # Dendritic cells
#     ["PPBP", "PF4"]                         # Megakaryocytes
# ]

# Flatten marker genes list
# marker_genes = [gene for sublist in marker_genes_list for gene in sublist]

# sc.pl.dotplot(
#     adata,
#     gene_marker_dict,
#     groupby=f"leiden_{chosen_res}",
#     dendrogram=True,
#     save="_marker_genes_dotplot11.png"
#     )

sc.tl.rank_genes_groups(
    adata,
    groupby=f"leiden_{chosen_res}",
    method="wilcoxon",
    n_genes=5,
    use_raw=True,
    key_added="rank_genes_groups"
    )

sc.pl.rank_genes_groups_dotplot(
    adata,
    # var_names=gene_marker_dict,
    n_genes=5,
    groupby=f"leiden_{chosen_res}",
    # dendrogram=True,
    save="_top5_genes.png"
)

# sc.tl.filter_rank_genes_groups(
#     adata,
#     groupby=f"leiden_{chosen_res}",
#     min_in_group_fraction=0.2, # only keep genes that are expressed in at least 25% of the cells in the cluster
#     max_out_group_fraction=0.4, # only keep genes that are expressed in less than 50% of the cells in other clusters
#     key="rank_genes_groups",
#     key_added="filtered_rank_genes_groups"
#     )

# sc.pl.rank_genes_groups_dotplot(
#     adata,
#     key="filtered_rank_genes_groups",
#     n_genes=5,
#     groupby=f"leiden_{chosen_res}",
#     standard_scale="var",
#     save="_filtered_top5_genes_dotplot11.png"
# )

# does this give dot plot or heatmap?
# sc.pl.rank_genes_groups(
#     adata,
#     n_genes=5,
#     # sharey=False,
#     save="_marker_genes11.png"
#     )

# Get top 5 genes across all clusters?
# top_genes = [] #flat
# names= adata.uns['rank_genes_groups']['names']
# print(names.dtype)
# print(names)

# # Remove duplicates
# top_genes = list(dict.fromkeys(top_genes))

# # Dotplot of top 5 marker genes per cluster
# sc.pl.dotplot(
#     adata,
#     top_genes,
#     groupby=f"leiden_{chosen_res}",
#     dendrogram=True,
#     save="_top5_genes_dotplot11.png"
# )

# Add after: adata = adata[~adata.obs.outlier & ~adata.obs.mt_outlier].copy()

# print("\n\nQC metrics by cluster (5 MADs filtering):")
# print("\nCluster 5 statistics:")
# cluster_5_mask = adata.obs[f"leiden_{chosen_res}"] == "5"
# print(f"  Cells in cluster 5: {cluster_5_mask.sum()}")
# print(f"  Mean total_counts: {adata.obs.loc[cluster_5_mask, 'total_counts'].mean():.1f}")
# print(f"  Median total_counts: {adata.obs.loc[cluster_5_mask, 'total_counts'].median():.1f}")
# print(f"  Mean n_genes_by_counts: {adata.obs.loc[cluster_5_mask, 'n_genes_by_counts'].mean():.1f}")
# print(f"  Mean pct_counts_mt: {adata.obs.loc[cluster_5_mask, 'pct_counts_mt'].mean():.1f}")

# print("\nOverall dataset statistics:")
# print(f"  Mean total_counts: {adata.obs['total_counts'].mean():.1f}")
# print(f"  Median total_counts: {adata.obs['total_counts'].median():.1f}")
# print(f"  Mean n_genes_by_counts: {adata.obs['n_genes_by_counts'].mean():.1f}")
# print(f"  Mean pct_counts_mt: {adata.obs['pct_counts_mt'].mean():.1f}")

# 12. Assign cell type labels to each cluster
# ---------------------------------------------------------------------
# print("\n\nQuestion 12: Assign cell type labels to each cluster")
# # Map leiden clusters to cell type labels
# # Based on resolution=0.50, we have 9 clusters.
# cluster_labels = {
#     "0": "CD4+ T cells", # high CD3D, no CD8A
#     "1": "CD8+ T cells", # small CD8A  and CD8B (highest CD8A of all clusters) and no CD4
#     "2": "CD4+ T cells", # high CD3D, no CD8A
#     "3": "FCGR3A+ Monocytes", # high FCGR3A and moderate MS4A7, low CD14. (highest MS4A7 of all clusters)
#     "4": "CD14+ Monocytes", # small/moderate CD14, with high LYZ and CST3. (highest CD14 of all clusters)
#     "5": "Dendritic cells", # high FCER1A and CST3, low CD14 
#     "6": "CD4+ T cells", # high CD3D, no CD8A,
#     "7": "B cells", # high MS4A1, CD79A, CD79B
#     "8": "NK cells" # high NKG7, GNLY and moderate KLRD1
# }

# adata.obs["cell_type"] = adata.obs[f"leiden_{chosen_res}"].map(cluster_labels).astype("category")

# sc.pl.dotplot(
#     adata,
#     gene_marker_dict,
#     groupby="cell_type",
#     dendrogram=True,
#     save="_marker_genes_dotplot_cell_type12.png"
#     )

# ax = sc.pl.heatmap(
#     adata,
#     gene_marker_dict,
#     groupby=f"leiden_{chosen_res}",
#     dendrogram=True,
#     save="_marker_genes_heatmap12.png",
#     cmap="viridis"
# )

# # 13. Plot UMAP with clusters colored by assigned cell type labels
# # ---------------------------------------------------------------------
# print("\n\nQuestion 13: Plot UMAP with clusters colored by assigned cell type labels")
# # after assigning cell type labels to each cluster, create a new column in the adata.obs dataframe that contains the cell type labels for each cell. 

# with plt.rc_context({"figure.figsize": (8,8)}):
#     sc.pl.umap(adata, 
#                color="cell_type", 
#                save="_cell_types13.png",
#                legend_loc="on data",
#                legend_fontsize=12,
#                legend_fontoutline=2,
#                alpha=0.5,
#                frameon=False,
#                title="UMAP colored by assigned cell type labels") 

# # Save the normalized and annotated AnnData object for future use
# ##############################################################################
# # adata.write(f"{OUTPUT_DIR}/adata.h5ad")

# # Close the log file
# log_file.close()

