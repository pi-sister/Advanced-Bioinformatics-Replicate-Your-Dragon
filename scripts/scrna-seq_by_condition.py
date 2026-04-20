# All mitochondrial QC was commented out, since the reference does not have mitochondrial genes
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
import warnings
from matplotlib.lines import Line2D
from scipy.stats import pointbiserialr
warnings.filterwarnings('ignore')
seed(5042026)

# User-defined parameters
PROJECT_DIR = "/work/TALC/mdsc519_2026w/students/jamie/Dragon"
CK_DIR = os.path.join(PROJECT_DIR, "data", "output", "cellranger", "SRR24952454")
POST_DIR = os.path.join(PROJECT_DIR, "data", "output", "cellranger", "SRR24952453")
CELLRANGER_OUTDIR = "/work/TALC/mdsc519_2026w/students/jamie/Dragon/data/output/cellranger"
# FIG_OUTPUTS_DIR = os.path.join(PROJECT_DIR, "output","MAD3PC20")

# Configure scanpy settings
sc.set_figure_params(figsize=(8, 8), dpi=100)
plt.rcParams["figure.figsize"] = (8, 8)
plt.rcParams["figure.dpi"] = 100

samples = {
    "post": "SRR24952453", "ck": "SRR24952454"
}

n_PCs = 20
chosen_res = 1.0
filtering_method = "IQR"
multiplier = 4.5

# Functions for QC and filtering
def get_expected_doublet_rate(n_cells):
    # Find the closest key in the lookup dictionary
    if n_cells <= 500:
        return 0.004
    elif n_cells <= 1000:
        return 0.008
    elif n_cells <= 2000:
        return 0.016
    elif n_cells <= 3000:
        return 0.023
    elif n_cells <= 4000:
        return 0.031
    elif n_cells <= 5000:
        return 0.039
    elif n_cells <= 6000:
        return 0.046
    elif n_cells <= 7000:
        return 0.054
    elif n_cells <= 8000:
        return 0.061
    elif n_cells <= 9000:
        return 0.069
    else:
        return 0.076
    

# Functions for QC and filtering
def per_sample_qc(adata, sample_name):
    print(f"\nQC metrics for sample {sample_name}:")
    print(f"  Number of cells: {adata.n_obs}")
    print(f"  Number of genes: {adata.n_vars}")
    # adata.var["mt"] = adata.var_names.str.startswith(("MT-","mt-"))
    sc.pp.calculate_qc_metrics(
        adata,
        log1p=True,
        inplace=True,
        percent_top=[20, 50, 100]
    )
    
    with plt.rc_context({"figure.figsize": (15, 6)}):
        sc.pl.violin(
            adata, 
            # ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
            ["total_counts", "n_genes_by_counts"],
            # log=True,
            jitter=0.4,
            multi_panel=True,
            stripplot=True,
            size=3,
            palette="Set2",
            save=f"_qc_{sample_name}.png"
        )
        plt.title("Total UMI Counts per Cell", fontsize=14)
        plt.ylabel("log10(Total UMI counts)", fontsize=20)
        plt.xlabel("", fontsize=20)
        
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save=f"_total_counts_vs_n_genes_{sample_name}.png")
    # sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save=f"_total_counts_vs_pct_counts_mt_{sample_name}.png")
    
    # print(f"Figures saved for sample {sample_name} QC metrics.")
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
        # IQR 3 is more stringent than IQR 1.5, as it allows for a wider range of values to be considered non-outliers. The choice of multiplier depends on the desired balance between retaining more cells and removing potential outliers. Using a multiplier of 3 may be appropriate if you want to be more conservative in filtering and retain more cells, while a multiplier of 1.5 may be suitable if you want to be more aggressive in filtering and remove more potential outliers.
        # for metric in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        for metric in ["total_counts", "n_genes_by_counts"]:
            lower_bound, upper_bound = iqr_bounds(filtered.obs[metric], multiplier)
            # print(f"  {metric}: [{lower_bound:.3f}, {upper_bound:.3f}]")
            filtered = filtered[
                (filtered.obs[metric] >= lower_bound) &
                (filtered.obs[metric] <= upper_bound)
            ].copy()
        # print(f"  Filtered number of cells: {filtered.n_obs}")
        return filtered
    # 
    elif isinstance(multiplier, (int, float)):
        # more systematic approach based on MAD (Median Absolute Deviation) to identify outliers in the QC metrics and filter out low-quality cells.
        for metric in ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_in_top_20_genes"]:
            M = filtered.obs[metric]
            median = np.median(M)
            mad = median_abs_deviation(M)
            lower = median - multiplier * mad
            upper = median + multiplier * mad
            print(f"  {metric}: [{lower:.3f}, {upper:.3f}]")

        filtered.obs["outlier"] = (
            is_outlier(filtered, "log1p_total_counts", multiplier)
            | is_outlier(filtered, "log1p_n_genes_by_counts", multiplier)
            | is_outlier(filtered, "pct_counts_in_top_20_genes", multiplier)
        )
        
    elif isinstance(multiplier, list):
        for pos, metric in enumerate(["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_in_top_20_genes"]):
            M = filtered.obs[metric]
            median = np.median(M)
            mad = median_abs_deviation(M)
            lower = median - multiplier[pos] * mad
            upper = median + multiplier[pos] * mad
            print(f"  {metric}: [{lower:.3f}, {upper:.3f}]")

        filtered.obs["outlier"] = (
            is_outlier(filtered, "log1p_total_counts", multiplier[0])
            | is_outlier(filtered, "log1p_n_genes_by_counts", multiplier[1])
            | is_outlier(filtered, "pct_counts_in_top_20_genes", multiplier[2])
        )
        
    print(f"  Outliers based on total counts, n_genes_by_counts, and pct_counts_in_top_20_genes: {filtered.obs.outlier.value_counts()}")
        
        # filtered.obs["mt_outlier"] = is_outlier(filtered, "pct_counts_mt", 3) | (
        #     filtered.obs["pct_counts_mt"] > 8
        # )                   
        # print(f"  Outliers based on pct_counts_mt: {filtered.obs.mt_outlier.value_counts()}")
        
        # filtered = filtered[~filtered.obs.outlier & ~filtered.obs.mt_outlier].copy()
    filtered = filtered[~filtered.obs.outlier].copy()
    print(f"  Filtered number of cells ({sample_name}): {filtered.n_obs}")
    return filtered
    
# scrublet to identify doublets - not sure if we have to do this since we already have the cellranger output, which should have already filtered out doublets. But could be worth trying to see if we can identify any additional doublets that were missed by cellranger.
def identify_doublets(adata, expected_doublet_rate=0.06, sample_name=None, umap=False):
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets
    adata.uns["scrublet"] = scrub
    print(f"Predicted doublets: {predicted_doublets.sum()} out of {adata.n_obs} cells")
    if sample_name is not None:
        print(f"DEBUG: Doublet score distribution for sample {sample_name}")
        # sc.pl.scrublet_score_distribution(adata, save=f"_doublet_score_distribution_{sample_name}.png")
        scrub.plot_histogram()
        plt.savefig(f"{fig_output_dir}/doublet_score_distribution_{sample_name}.png", bbox_inches='tight', dpi=150)
        plt.close()
    if umap:
        sc.pl.umap(adata, color=['doublet_score', 'predicted_doublet'], save=f"_doublet_{sample_name}.png")
    # TODO: finish the inspection and choice of doublet threshold
#     adata.obs["possible_doublets"] = adata.obs["predicted_doublet"].copy()
    
#     # stricter threshold for doublets
#     inspect = adata.obs["possible_doublets"].values
#     inspect_scores = adata.obs.loc[inspect, "doublet_score"]
    
#     strict_cut = np.quantile(inspect_scores, 0.5)
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    return adata
    
    
def pct_expressed_by_group(adata, groupby, layer=None):
    """ Helper function to calculate the in and out frequency by cluster

    Args:
        adata (AnnData): expression object
        groupby (string): variable name to group by
        layer (Bool, optional): Layer of AnnData with X. Defaults to None.

    Returns:
        dict: with percents in and out across all genes for each cluster
    """
    X = adata.layers[layer] if layer is not None else adata.X
    if sparse.issparse(X):
        X = X.tocsr()

    # groups is the cluster labels for each cell
    groups = adata.obs[groupby].astype(str)
    result = {}

    # for each cluster
    for g in sorted(groups.unique()):
        # index of cells in the cluster and outside the cluster
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


#########################################################################################
# MAIN ANALYSES
#########################################################################

# Redirect all output to a log file
# log_file = open(f"{OUTPUT_DIR}/output.log", "w")
# sys.stdout = log_file

adatas = []
adata_markers_list = []
for sample, srr_id in samples.items():
    settings_string = f"PCs{n_PCs}_res{chosen_res}_{filtering_method}{multiplier}"
    fig_output_dir = os.path.join(PROJECT_DIR, "final_output", sample)
    os.makedirs(fig_output_dir, exist_ok=True)
    sc.settings.figdir = os.path.join(fig_output_dir)
    h5_path = os.path.join(CELLRANGER_OUTDIR, sample, "outs", "filtered_feature_bc_matrix.h5")
    adata = sc.read_10x_h5(h5_path)
    
    # Use the feature names from the H5 file (usually the second column contains gene symbols)
    # Check what columns exist and use the appropriate one
    if "gene_symbols" in adata.var.columns:
        adata.var_names = adata.var["gene_symbols"].values
    elif "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].values
    else:
        # If neither exists, keep the default (usually gene_ids)
        print(f"  Available var columns: {adata.var.columns.tolist()}")
    
    adata.var_names_make_unique() # make gene names unique by adding a suffix to duplicate names, which is important for downstream analyses that require unique gene identifiers. This step ensures that each gene can be uniquely identified and prevents issues that may arise from having duplicate gene names in the dataset.
    adata.obs["sample_id"] = srr_id
    adata.obs["condition"] = sample

    # for each adata, comput qc, remove outliers. Then concatenate the adatas together for downstream analyses. 
    adata = per_sample_qc(adata, sample)
    
    print(f"\nCells before filtering: {adata.n_obs}")
    # Use 3 for more stringent filtering, 1.5 for less stringent filtering. The paper does not specify the exact multiplier used for the IQR method, so we can experiment with different values to see how it affects the number of cells retained and the downstream analyses.
    # adata = per_sample_filtering(adata, sample, method="IQR", multiplier=3) 
    # MAD filtering. Use multiplier of 3 for less deviation from the median, 5 for more deviation. 
    adata = per_sample_filtering(adata, sample, method=filtering_method, multiplier=multiplier)
    adata.obs["filtered"] = True
    
    expected_doublet_rate = get_expected_doublet_rate(adata.n_obs)
    adata = identify_doublets(adata, expected_doublet_rate=expected_doublet_rate, sample_name=sample)
    
    adatas.append(adata)

    # adata = sc.concat(adatas, join="outer", label="sample_id", keys=list(samples.keys()))
    adata.obs_names_make_unique()

    # Normalize the data
    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)


    # Identify highly variable genes
    sc.pp.highly_variable_genes(
        adata, # use the filtered adata?  to identify highly variable genes, as the normalization step does not change the variance of the genes, so we can use the original data to identify highly variable genes. (Get a warning when using normalized data) - flavor seurat expects normalized data, seurat_v3 expects raw counts data
        flavor="seurat",
        batch_key="condition",
        # span=0.5, 
        # n_top_genes=200,
        # inplace=False
        # min_mean=0.05, 
        # max_mean=1, 
        # min_disp=0.8,
        )
    sc.pl.highly_variable_genes(adata, save=".png")
    print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")

    adata_markers = adata.copy()
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)


    #########################################################################################
    # Dimensionality Reduction and Clustering
    #########################################################################################
    sc.pp.pca(adata, n_comps=50, svd_solver="arpack")
    sc.pl.pca_variance_ratio(adata, log=True, save=".png")
    sc.pl.pca_overview(adata, save=".png")

    # Elbow plot to show variance explained by each principal component (first 40 PCs)
    variance = adata.uns['pca']['variance'][:40]
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(variance)+1), variance, 'bo-', linewidth=2, markersize=6)
    plt.xlabel('Principal Component', fontsize=20)
    plt.ylabel('Variance Explained', fontsize=20)
    plt.title('Elbow Plot', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{fig_output_dir}/pca_elbow_plot.png", bbox_inches='tight')
    plt.close()

    n_PCs = n_PCs # choose the number of PCs to use for downstream analyses

    # Investigate PC - correlation with condition and QC metrics

    pcs = adata.obsm["X_pca"]
    cond = (adata.obs["condition"] == "post").astype(int).values

    print("Correlation of PCs with condition (point biserial correlation):")
    for i in range(10):
        r, p = pointbiserialr(cond, pcs[:, i])
        print(f"PC{i+1}: r={r:.3f}, p={p:.3e}")

    for col in ["total_counts", "n_genes_by_counts", "pct_counts_in_top_20_genes"]:
        print(f"\nCorrelation with {col}")
        vals = adata.obs[col].values
        for i in range(10):
            r = np.corrcoef(vals, pcs[:, i])[0, 1]
            print(f"PC{i+1}: r={r:.3f}")
            
    # Construct kNN graph on pca and perform clustering using the Leiden algorithm
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_PCs)
    resolutions = [0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
    res_table = []
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", n_iterations=2, flavor="igraph", directed=False)
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

    # Choose resolution and compute UMAP
    chosen_res = chosen_res 
    groups = adata.obs[f"leiden_{chosen_res}"]
    print(f"Resolution: {chosen_res}, Clusters: {groups.nunique()}")

    sc.tl.umap(adata)
    # colour by cluster
    # sc.pl.umap(adata, 
    #         color=f"leiden_1.0", 
    #     #    add_outline=True,
    #         legend_loc="on data",
    #         legend_fontsize=20,
    #         alpha=0.5,
    #     #    legend_fontoutline=2,
    #         frameon=False,
    #         title=f"UMAP colored by Leiden clusters (resolution=1.0)",
    #         save="_clustering9.png")
    # Plot UMAP colored by different metrics
    sc.pl.umap(adata, color="total_counts", frameon=False, save="_total_counts.png")
    # sc.pl.umap(adata, color="pct_counts_mt", frameon=False, save="_pct_counts_mt.png")
    sc.pl.umap(adata, color=f"leiden_{chosen_res}", alpha=0.5, frameon=False, legend_loc="on data", legend_fontsize=20, title=f"UMAP colored by Leiden clusters (resolution={chosen_res})", save=".png")

    # give cell cluster labels to the adata_markers object for marker discovery
    adata_markers.obs[f"leiden_{chosen_res}"] = adata.obs[f"leiden_{chosen_res}"].values


    # Cluster Annotation
    #########################################################################################

    # Identify marker genes for each cluster
    sc.tl.dendrogram(adata_markers, groupby=f"leiden_{chosen_res}")

    sc.tl.rank_genes_groups(
        adata_markers,
        groupby=f"leiden_{chosen_res}",
        method="wilcoxon",
        use_raw=False,
        pts=False, # manually calculate pct_in and pct_out later, as the default implementation does not work with sparse matrices and is very slow for large datasets.
        key_added="cluster_markers" # store the results in adata.uns["cluster_markers"]
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
        
        # genes filtered by position
        gene_to_idx = {gene: idx for idx, gene in enumerate(adata_markers.var_names)}
        df["pct_in"] = df["names"].map(lambda gene: pct_info[cluster]["pct_in"][gene_to_idx[gene]])
        df["pct_out"] = df["names"].map(lambda gene: pct_info[cluster]["pct_out"][gene_to_idx[gene]])
        df["cluster"] = cluster
        df["gene_idx"] = df["names"].map(gene_to_idx)
        all_markers.append(df)
        
    markers_df = pd.concat(all_markers, ignore_index=True)
    markers_df.to_csv(os.path.join(fig_output_dir, f"markers_with_pct_info.csv"), index=False)

    # Thresholds from the paper!
    filtered_markers_df = markers_df[
        (markers_df["pct_in"] > 0.4) &
        (markers_df["pct_out"] < 0.2) &
        (markers_df["pvals_adj"] <= 0.01) &
        (markers_df["logfoldchanges"] >= 1.0)
    ].copy()

    # Rank and keep top 3 per cluster
    top3_markers = (
        filtered_markers_df.sort_values(["cluster", "logfoldchanges"], ascending=[True, False])
        .groupby('cluster', group_keys=False).head(3).reset_index(drop=True)
    )

    filtered_markers_df.to_csv(os.path.join(fig_output_dir, f"filtered_markers_with_pct_info.csv"), index=False)
    top3_markers.to_csv(os.path.join(fig_output_dir, f"top3_markers.csv"), index=False)

    print(f"Top 3 markers per cluster: \n{top3_markers[['cluster', 'names', 'logfoldchanges', 'pct_in', 'pct_out']]}")


    # Plot expression of top 3 markers per cluster
    # ---------------------------------------------------------------------
    markers_list = top3_markers["names"].tolist()

    # axes = sc.pl.dotplot(
    #     adata_markers, 
    #     var_names=markers_list, 
    #     groupby=f"leiden_{chosen_res}", 
    #     standard_scale="var",
    #     # use_raw=False, 
    #     dendrogram=False,
    #     save="_top3_markers.png"
    # )

    cluster_col = f"leiden_{chosen_res}"

    # numeric cluster order
    cluster_order = sorted(adata_markers.obs[cluster_col].astype(int).unique())
    cluster_order_str = [str(x) for x in cluster_order]

    # enforce plotting order
    adata_markers.obs[cluster_col] = pd.Categorical(
        adata_markers.obs[cluster_col].astype(str),
        categories=cluster_order_str,
        ordered=True
    )

    # sort top markers by numeric cluster, not string cluster
    top3_markers = top3_markers.copy()
    top3_markers["cluster"] = top3_markers["cluster"].astype(str)
    top3_markers["cluster_num"] = top3_markers["cluster"].astype(int)

    top3_markers = top3_markers.sort_values(
        ["cluster_num", "logfoldchanges"],
        ascending=[True, False]
    )

    # grouped var_names
    top3_dict = (
        top3_markers.groupby("cluster")["names"]
        .apply(list)
        .reindex(cluster_order_str)
        .to_dict()
    )

    # Remove NaN values from the dictionary to avoid "float object is not iterable" error
    top3_dict = {k: v for k, v in top3_dict.items() if isinstance(v, list) and len(v) > 0}

    sc.pl.dotplot(
        adata_markers,
        var_names=top3_dict,
        groupby=cluster_col,
        standard_scale="var",
        use_raw=False,
        save="_top3_markers_.png"
    )

    print(top3_dict)


    # Plot expression of top 3 markers per cluster on UMAP
    top3_genes = top3_markers["names"].unique().tolist()

    # Create a grid of UMAP plots for each gene
    n_genes = len(top3_genes)
    n_cols = 5
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4*n_rows), dpi=100)
    axes = axes.flatten() if n_genes > 1 else np.array([axes])

    for idx, gene in enumerate(top3_genes):
        try:
            sc.pl.umap(
                adata_markers,
                color=gene,
                frameon=True,
                title=f"{gene}",
                ax=axes[idx],
                show=False,
                cmap="viridis"
            )
        except Exception as e:
            print(f"Error plotting gene {gene}: {e}")
            axes[idx].text(0.5, 0.5, f"Error plotting {gene}", ha='center', va='center')

    # Hide unused subplots
    for idx in range(n_genes, len(axes)):
        axes[idx].axis('off')

    plt.tight_layout()
    plt.savefig(f"{fig_output_dir}/umap_top3_markers.png", bbox_inches='tight', dpi=150)
    plt.close()

    marker_panels = {
        "Mesocarp": [
            "HU08G02237", "HU11G01114", "HU09G01127", "HU06G01840",
            "HU10G01409", "HU08G00805", "HU02G02701", "HU01G00649"
        ],
        "Exocarp": [
            "HU07G01483", "HU03G02606", "HU01G01880",
            "HU09G00039", "HU08G01941", "HU07G01714",
            "HU03G00170", "HU06G01661", "HU08G01266"
        ],
        "Endocarp": [
            "HU03G02180", "HU03G02824", "HU07G01763",
            "HU01G02711", "HU08G01810", "HU06G02555"
        ],
        "Endocarp fiber": [
            "HU05G00062", "HU05G00061", "HU10G00161", "HU07G02077"
        ],
        "Vascular bundle": [
            "HU05G00061", "HU05G01893", "HU01G01937", "HU10G00163"
        ],
    }

    highlighted_genes = ["HU08G01266", "HU06G02555", "HU07G02077", "HU08G02237", "HU10G00163"]

    # sc.tl.dendrogram(adata_markers, groupby=f"leiden_{chosen_res}")
    axes = sc.pl.dotplot(
        adata_markers, 
        marker_panels, 
        groupby=f"leiden_{chosen_res}", 
        # swap_axes=True,
        # standard_scale="var",
        use_raw=False,
        show=False,
        cmap="Greens"
    )
    print("Dotplot axes:", axes.keys())
    ax = axes["mainplot_ax"]
    for label in ax.get_xticklabels():
        gene_name = label.get_text()
        if gene_name in highlighted_genes:
            label.set_color("#c60d61")
            label.set_fontweight('bold')
            
    # Add a legend for the highlighted genes
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='ST Identified', markerfacecolor='#c60d61', markersize=10)]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1))
    plt.tight_layout()
    plt.savefig(f"{fig_output_dir}/dotplot_highlighted.png", bbox_inches='tight', dpi=150)
    plt.close()

    # Assign cell types
    # ---------------------------------------------------------------------
    ck_cluster_to_celltype = {
        "0": "EX",
        "1": "EX",
        "2": "EX",
        "3": "EX",
        "4": "EX",
        "5": "EN-LIKE",
        "6": "EX",
        "7": "EX",
        "8": "EX-LIKE",
        "9": "EX-EN-MIX",
        "10": "EN",
        "11": "VAS",
        "12": "EX",
        "13": "EN",
    }
    
    post_cluster_to_celltype = {
        "0": "ME",
        "1": "ME",
        "2": "ME",
        "3": "ME",
        "4": "ME",
        "5": "ME",
        "6": "ME",
        "7": "ME",
        "8": "VAS",
        "9": "ME",
        "10": "ME",
        "11": "EX-LIKE",
        "12": "ME",
        "13": "EN",
        "14": "EN",
        "15": "EN-LIKE",
    }

    if sample == "ck":
        cluster_to_celltype = ck_cluster_to_celltype
    else:
        cluster_to_celltype = post_cluster_to_celltype

    adata_markers.obs["cell_type"] = adata_markers.obs[f"leiden_{chosen_res}"].map(cluster_to_celltype)

    # Copy UMAP embedding from adata to adata_markers if it exists
    if 'X_umap' in adata.obsm:
        adata_markers.obsm['X_umap'] = adata.obsm['X_umap']

    sc.pl.umap(adata_markers, alpha=0.5, color="cell_type", frameon=False, legend_loc="on data", legend_fontsize=20, title=f"UMAP colored by cell type", save="_cell_type.png")

    # Append adata_markers with cell type annotations to the list for downstream analysis
    adata_markers_list.append(adata_markers)

    # Compare post and ck conditions
    # -------------------------------------------------------------------
    # sc.pl.umap(adata_markers, color="condition", frameon=False, legend_loc="on data", legend_fontsize=20, ncols=2, save="_condition.png")

    # post_ad = adata_markers[adata_markers.obs["condition"] == "post"].copy()
    # ck_ad = adata_markers[adata_markers.obs["condition"] == "ck"].copy()

    # #scanpy example on the "customizing scanpy plots" page
    # ncols=2
    # nrows=1
    # figsize=4
    # wspace=0.5
    # fig, axs = plt.subplots(
    #     nrows=nrows,
    #     ncols=ncols,
    #     figsize=(ncols*figsize+figsize+wspace+(ncols-1), nrows*figsize)
    # )

    # plt.subplots_adjust(wspace=wspace)
    # sc.pl.umap(
    #     post_ad, 
    #     color=f"leiden_{chosen_res}", 
    #     frameon=False, 
    #     legend_loc="on data", 
    #     legend_fontsize=20, 
    #     title=f"Post condition",
    #     ax=axs[0],
    #     show=False
    # )

    # sc.pl.umap(
    #     ck_ad, 
    #     color=f"leiden_{chosen_res}", 
    #     frameon=False, 
    #     legend_loc="on data", 
    #     legend_fontsize=20, 
    #     title=f"CK condition",
    #     ax=axs[1],
    #     show=False
    # )

    # plt.savefig(f"{fig_output_dir}/umap_comparison.png", bbox_inches='tight', dpi=150)
    plt.close()
    
    
fig_output_dir = os.path.join(PROJECT_DIR, "final_output")

marker_panels = {
    "ME": ["HU08G02237", "HU11G01114", "HU09G01127", "HU06G01840",
           "HU10G01409", "HU08G00805", "HU02G02701", "HU01G00649"],
    "EX": ["HU07G01483", "HU03G02606", "HU01G01880",
           "HU09G00039", "HU08G01941", "HU07G01714",
           "HU03G00170", "HU06G01661", "HU08G01266"],
    "EN": ["HU03G02180", "HU03G02824", "HU07G01763",
           "HU01G02711", "HU08G01810", "HU06G02555"],
    "ENF": ["HU05G00062", "HU05G00061", "HU10G00161", "HU07G02077"],
    "VB": ["HU05G00061", "HU05G01893", "HU01G01937", "HU10G00163"],
}

def cluster_panel_scores(adata, cluster_col):
    out = []
    for cl in sorted(adata.obs[cluster_col].astype(str).unique(), key=int):
        sub = adata[adata.obs[cluster_col].astype(str) == cl].copy()
        row = {"cluster": cl, "n_cells": sub.n_obs}
        for panel_name, genes in marker_panels.items():
            genes_present = [g for g in genes if g in sub.var_names]
            if len(genes_present) == 0:
                row[panel_name] = np.nan
            else:
                X = sub[:, genes_present].X
                if hasattr(X, "toarray"):
                    X = X.toarray()
                row[panel_name] = X.mean()
        out.append(row)
    df = pd.DataFrame(out)
    df["best_panel"] = df[["ME", "EX", "EN", "ENF", "VB"]].idxmax(axis=1)
    return df

# Example for one sample object
ck_scores = cluster_panel_scores(adata_markers_list[1], "leiden_1.0")
post_scores = cluster_panel_scores(adata_markers_list[0], "leiden_1.0")

print(ck_scores)
print(post_scores)

ck_scores.to_csv(os.path.join(fig_output_dir, f"ck_cluster_panel_scores.csv"), index=False)
post_scores.to_csv(os.path.join(fig_output_dir, f"post_cluster_panel_scores.csv"), index=False)


all_adatas = sc.concat(adata_markers_list, join="outer", label="sample_id", keys=list(samples.keys()))
# all_datas["ck"] = adata_markers_list[1]
# all_datas["post"] = adata_markers_list[0]
all_adatas.obs_names_make_unique()


de_results = {}
for ct in all_adatas.obs["cell_type"].dropna().unique():
    sub = all_adatas[all_adatas.obs["cell_type"] == ct].copy()
    sub.obs["condition"] = sub.obs["condition"].astype(str) 
    print(f"\nProcessing cell type: {ct}, cells: {sub.n_obs}")

    cond_counts = sub.obs["condition"].value_counts()
    print(f"\n=== {ct} ===")
    print(cond_counts)

    if "ck" not in cond_counts.index or "post" not in cond_counts.index:
        print("Skipping: missing ck or post")
        continue

    if cond_counts["ck"] < 3 or cond_counts["post"] < 3:
        print("Skipping: too few cells in one condition")
        continue

    ct_safe = str(ct).replace(" ", "_").replace("/", "_")
    key_name = f"de_{ct_safe}"

    # DE: senescent vs mature
    sc.tl.rank_genes_groups(
        sub,
        groupby="condition",
        groups=["post"],
        reference="ck",
        method="wilcoxon",
        pts=True,
        key_added=key_name
    )

    de_df = sc.get.rank_genes_groups_df(sub, group="post", key=key_name)
    de_df["cell_type"] = ct
    de_results[ct] = de_df

    # save table
    de_df.to_csv(f"{fig_output_dir}/de_{ct_safe}_post_vs_ck.csv", index=False)

    # --------------------
    # volcano plot
    # --------------------
    df = de_df.copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchanges", "pvals_adj"])
    df["pvals_adj"] = df["pvals_adj"].replace(0, 1e-300)
    df["neglog10_padj"] = -np.log10(df["pvals_adj"])

    lfc_thresh = 0.25
    p_thresh = 0.05

    df["status"] = "not_sig"
    df.loc[(df["logfoldchanges"] >= lfc_thresh) & (df["pvals_adj"] <= p_thresh), "status"] = "up"
    df.loc[(df["logfoldchanges"] <= -lfc_thresh) & (df["pvals_adj"] <= p_thresh), "status"] = "down"

    plt.figure(figsize=(7, 6))
    for status in ["not_sig", "up", "down"]:
        part = df[df["status"] == status]
        plt.scatter(part["logfoldchanges"], part["neglog10_padj"], s=18, alpha=0.7, label=status)

    plt.axvline(lfc_thresh, linestyle="--", linewidth=1)
    plt.axvline(-lfc_thresh, linestyle="--", linewidth=1)
    plt.axhline(-np.log10(p_thresh), linestyle="--", linewidth=1)

    sig = df[df["status"] != "not_sig"].copy()
    label_df = pd.concat([
        sig[sig["status"] == "up"].nlargest(6, "neglog10_padj"),
        sig[sig["status"] == "down"].nlargest(6, "neglog10_padj")
    ]).drop_duplicates()

    for _, row in label_df.iterrows():
        plt.text(row["logfoldchanges"], row["neglog10_padj"], row["names"], fontsize=8)

    plt.xlabel("log2 fold change")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(f"{ct}: post vs ck")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"{fig_output_dir}/volcano_{ct_safe}.png", bbox_inches='tight', dpi=150)
    plt.show()

    # --------------------
    # heatmap of top DE genes
    # --------------------
    up = df[(df["logfoldchanges"] >= lfc_thresh) & (df["pvals_adj"] <= p_thresh)] \
        .sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, False]) \
        .head(6)

    down = df[(df["logfoldchanges"] <= -lfc_thresh) & (df["pvals_adj"] <= p_thresh)] \
        .sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, True]) \
        .head(6)

    top_genes = up["names"].tolist() + down["names"].tolist()
    top_genes = [g for g in top_genes if g in sub.var_names]

    if len(top_genes) > 0:
        sc.pl.heatmap(
            sub,
            var_names=top_genes,
            groupby="condition",
            swap_axes=True,
            standard_scale="var",
            show_gene_labels=True,
            save=f"{ct_safe}.png"
        )

        # optional violin plots for a few top genes
        violin_genes = top_genes[:5]
        sc.pl.violin(
            sub,
            keys=violin_genes,
            groupby="condition",
            stripplot=False,
            jitter=0.1,
            multi_panel=True,
            save=f"{ct_safe}.png"
        )
    else:
        print("No significant genes found for heatmap/violin with current thresholds.")
    
    