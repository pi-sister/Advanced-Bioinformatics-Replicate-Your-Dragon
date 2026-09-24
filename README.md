# **Advanced Bioinformatics Final Project**

This was the final project completed for MDSC 519: Advanced Bioinformatics. It is a replication of single cell and spatial analytics of RNA from the dragonfruit.


## Background

 - Li et al. (1) uses scRNA-seq and ST to examine the changes in transcription in different cell types during senescence. 
   
 - Fruit senescence creates lost revenue and food waste. Understanding this process could help identify targets for improving shelf life. 

-  The original paper sought to identify cell types, expression patterns and spatial patterns associated with senescence.  
  - This study aimed to determine whether similar cluster structure, marker genes, and
   senescence-related patterns can be found with modified methods and the information provided.



1. Li X, Li B, Gu S, Pang X, Mason P, Yuan J, et al. Single-cell and spatial RNA sequencing reveal the spatiotemporal trajectories of fruit senescence. Nat Commun. 2024 Apr 10;15(1):3108.

## Method
  

This project analyzes single-cell RNA-seq and spatial transcriptomics (Visium) data from dragon fruit (pitaya) pericarp tissue, comparing control (**ck**) and treated (**post**) conditions using two public SRA BioProjects (`PRJNA974579`, `PRJNA1002459`). Raw sequencing data were retrieved with **NCBI Entrez Direct** (`esearch`/`efetch`) and the **SRA Toolkit** (`prefetch`, `fasterq-dump`), then compressed and renamed to 10x conventions with `pigz`. A custom reference index was built for the `PitayaGenomic` assembly using **Cell Ranger** (`mkref`, v8.0.1), and reads were aligned/counted per sample with **Cell Ranger count** (scRNA-seq) and **Space Ranger count** (v4.0.1, Visium spatial sample).

  

Downstream scRNA-seq analysis used **Scanpy**, following a standard preprocessing-to-clustering workflow: QC metrics (`calculate_qc_metrics`), outlier filtering (IQR- and MAD-based methods, mitochondrial filtering omitted since the reference lacks annotated MT genes), doublet detection with **Scrublet**, normalization (`normalize_total` + `log1p`), highly-variable-gene selection, scaling, PCA (with point-biserial/Pearson correlation checks against condition and QC covariates), a kNN graph, and **Leiden** clustering swept across resolutions (final resolution 1.0) with UMAP embedding. Cluster marker genes were identified via **Wilcoxon rank-sum tests** (`rank_genes_groups`), filtered by percent-expressed and fold-change thresholds, and used together with curated marker panels to manually annotate clusters into cell types (Mesocarp, Exocarp, Endocarp, Endocarp fiber, Vascular bundle). Differential expression between post and ck conditions was computed per cell type (Wilcoxon), visualized as volcano plots, heatmaps, and violin plots.

  

A parallel spatial transcriptomics script uses **Squidpy** for normalization, HVG selection, PCA, neighbors, **Louvain** clustering, and UMAP, and was used to cross-check pericarp marker genes against the scRNA-seq clusters (though this script currently points at Squidpy's demo dataset rather than the project's own Visium output). Two planned integration/atlas scripts (`integrate.py`, `analyze_atlas.py`) are present but not yet implemented. The full pipeline is orchestrated on an HPC cluster via **SLURM** batch scripts using conda-managed environments.


## Conclusions

The replication partially reproduced the original study. Broad cell-type marker patterns and the mature to senescent shift were recovered, supporting the idea that the original findings reflect real biological structure. The discrepancy in the clustering split suggests the missing details and data limited reproducibility.
