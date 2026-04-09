
import os
from random import seed
import scanpy as sc
import scrublet as scr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
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
    return adata
    
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
    
    
def pct_expressed_by_group(adata, groupby, layer=None):
    X = adata.layers[layer] if layer is not None else adata.X
    if sparse.issparse(X):
        X = X.tocsr()

    groups = adata.obs[groupby].astype(str)
    result = {}

    for g in sorted(groups.unique()):
        idx_in = np.where(groups.values == g)[0]
        idx_out = np.where(groups.values != g)[0]

        X_in = X[idx_in]
        X_out = X[idx_out]

        if sparse.issparse(X):
            pct_in = np.asarray((X_in > 0).mean(axis=0)).ravel()
            pct_out = np.asarray((X_out > 0).mean(axis=0)).ravel()
        else:
            pct_in = (X_in > 0).mean(axis=0)
            pct_out = (X_out > 0).mean(axis=0)

        result[g] = {
            "pct_in": pct_in,
            "pct_out": pct_out
        }

    return result


    
# Main analysis workflow
#########################################################################

# Redirect all output to a log file
# log_file = open(f"{OUTPUT_DIR}/output.log", "w")
# sys.stdout = log_file

adatas = []
for sample, srr_id in samples.items():
    h5_path = os.path.join(CELLRANGER_OUTDIR, sample, "outs", "filtered_feature_bc_matrix.h5")
    adata = sc.read_10x_h5(h5_path, var_names="gene_symbols")
    adata.var_names_make_unique() # make gene names unique by adding a suffix to duplicate names, which is important for downstream analyses that require unique gene identifiers. This step ensures that each gene can be uniquely identified and prevents issues that may arise from having duplicate gene names in the dataset.
    adata.obs["sample_id"] = srr_id
    adata.obs["condition"] = sample

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
adata_markers = adata.copy()
sc.pp.scale(adata, max_value=10) # scale the data to unit variance and zero mean, with a maximum value of 10 to prevent extreme values from dominating the analysis. This step is important because it ensures that all genes are on the same scale and that highly expressed genes do not dominate the downstream analyses such as PCA and clustering. Scaling makes analyses more sensitive to subtle differences in gene expression across cells.


################################################################################
# Dimensionality Reduction and Clustering
################################################################################
# Run PCA and visualize the variance explained by each principal component
# ---------------------------------------------------------------------
print("\n\nRun PCA and visualize the variance explained by each principal component")
sc.pp.pca(adata, n_comps=50, svd_solver="arpack") 
# sc.tl.pca(adata, svd_solver="arpack") This version is deprecated? use sc.pp.pca instead
sc.pl.pca_variance_ratio(adata, log=True, save=".png")

sc.pl.pca_overview(adata, save=".png")

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
print("\n\nChoose resolution")
chosen_res = 0.50
print(f"\nChosen resolution: {chosen_res}")
# print the number of clusters and their sizes for the chosen resolution
print(f"Number of clusters: {adata.obs[f"leiden_{chosen_res}"].value_counts()}")

# UMAP embedding and plot cells
# ---------------------------------------------------------------------
#https://scanpy.readthedocs.io/en/stable/tutorials/plotting/core.html
print("\n\nUMAP embedding and plot")
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
# Cluster Annotation TODO: marker genes needed
################################################################################
# followed: https://scanpy.readthedocs.io/en/stable/tutorials/plotting/core.html
# Identify marker genes for each cluster
# ---------------------------------------------------------------------
print("\n\nIdentify marker genes for each cluster")
sc.tl.dendrogram(adata, groupby=f"leiden_{chosen_res}")

sc.tl.rank_genes_groups(
    adata_markers,
    groupby=f"leiden_{chosen_res}",
    method="wilcoxon",
    use_raw=False,
    pts=False,
    key_added="cluster_markers"
    )

pct_info = pct_expressed_by_group(adata_markers, groupby=f"leiden_{chosen_res}")


# convert into one dataframe
all_markers = []
clusters = adata_markers.obs[f"leiden_{chosen_res}"].astype(str).unique().tolist()

for cluster in sorted(clusters):
    df = sc.get.rank_genes_groups_df(
        adata_markers, 
        group=cluster, 
        key="cluster_markers"
    ).copy()
    
    gene_to_idx = {gene: idx for idx, gene in enumerate(adata_markers.var_names)}
    df["pct_in"] = df["names"].map(lambda gene: pct_info[cluster]["pct_in"][gene_to_idx[gene]])
    df["pct_out"] = df["names"].map(lambda gene: pct_info[cluster]["pct_out"][gene_to_idx[gene]])
    df["cluster"] = cluster
    df["gene_idx"] = df["names"].map(gene_to_idx)
    all_markers.append(df)
    
markers_df = pd.concat(all_markers, ignore_index=True)
markers_df.to_csv(os.path.join(FIG_OUTPUTS_DIR, f"markers_with_pct_info.csv"), index=False)

# Thresholds from the paper!
filtered_markers_df = markers_df[
    (markers_df["pct_in"] > 0.4) &
    (markers_df["pct_out"] < 0.2) &
    (markers_df["pvals_adj"] <= 0.01) &
    (markers_df["logfoldchanges"] >= 1.0)
].copy()

# Rank and keep top 3 per cluster
top3_markers = (
    filtered_markers.sort_values(["cluster", "logfoldchanges"], ascending=[True, False]).head(3).reset_index(drop=True)
)

filtered_markers_df.to_csv(os.path.join(FIG_OUTPUTS_DIR, f"filtered_markers_with_pct_info.csv"), index=False)
top3_markers.to_csv(os.path.join(FIG_OUTPUTS_DIR, f"top3_markers.csv"), index=False)

print(f"Top 3 markers per cluster: \n{top3_markers[['cluster', 'names', 'logfoldchanges', 'pct_in', 'pct_out']]}")


# Plot expression of top 3 markers per cluster
# ---------------------------------------------------------------------
markers_list = top3_markers["names"].tolist()

sc.pl.dotplot(
    adata_markers, 
    var_names=markers_list, 
    groupby=f"leiden_{chosen_res}", 
    standard_scale="var",
    dendogram=False,
    use_raw=False, 
    save="_top3_markers.png"
)

top3_dict = (
    top3_markers.groupby("cluster")["names"]
    .apply(list)
    .to_dict()
)

print(top3_dict)