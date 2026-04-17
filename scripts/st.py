import os

import squidpy as sq
import matplotlib.pyplot as plt
from random import seed

seed(5042026)

PROJECT_DIR = "/work/TALC/mdsc519_2026w/students/jamie/Dragon"
ST_ANNDATA_DIR = os.path.join(PROJECT_DIR, "data", "output", "cellranger", "visium")
FIG_OUTPUTS_DIR = os.path.join(PROJECT_DIR, "output", "st")
N_TOP_GENES = 2000


adata = sq.datasets.visium_hne_adata()
adata.obs_names_make_unique()
adata.var_names_make_unique()

sq.settings.figdir = os.path.join(FIG_OUTPUTS_DIR)
sq.set_figure_params(figsize=(8, 8), dpi=100)
plt.rcParams["figure.figsize"] = (8, 8)
plt.rcParams["figure.dpi"] = 100

# tissue-spot filtering

# SCTransform normalization and dimensionality reduction
sq.pp.normalize_total(adata, target_sum=1e4)
sq.pp.log1p(adata)
sq.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=N_TOP_GENES)
adata = adata[:, adata.var.highly_variable].copy()

sq.tl.pca(adata, n_comps=50)
sq.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sq.tl.louvain(adata, resolution=0.5)
sq.tl.umap(adata)
sq.pl.umap(adata, color=["louvain"], save="_clusters.png")