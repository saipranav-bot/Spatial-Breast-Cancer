# Spatial Breast Cancer Transcriptomics Analysis

## Overview

A reproducible spatial transcriptomics analysis pipeline for breast cancer Visium data using Seurat.

This project performs quality control, normalization, dimensionality reduction, clustering, marker gene identification, differential expression analysis, functional enrichment analysis, spatial variable gene detection, and spatial annotation.

---

## Analysis Workflow

### 1. Data Loading

* Import 10x Genomics Visium spatial transcriptomics data
* Create Seurat object
* Prepare spatial expression matrix

### 2. Quality Control

* Feature count filtering
* UMI count filtering
* Mitochondrial content analysis
* QC visualization

### 3. Normalization

* SCTransform normalization
* Variance stabilization
* Feature selection

### 4. Dimensionality Reduction

* PCA analysis
* Elbow plot analysis
* Principal component selection

### 5. Spatial Clustering

* Graph-based clustering
* Spatial domain identification
* UMAP visualization

### 6. Marker Gene Analysis

* Cluster marker identification
* Heatmap visualization
* Dotplot visualization

### 7. Cell Type Annotation

* Marker-based annotation
* Spatial validation
* Annotated spatial maps

### 8. Differential Expression Analysis

* Cluster-specific differential expression analysis
* Identification of differentially expressed genes

### 9. Functional Enrichment Analysis

Performed using:

* Gene Ontology (GO)
* KEGG pathway enrichment

### 10. Spatial Variable Gene Analysis

* Identification of genes showing spatial expression patterns

---

# Tools and Software

## Programming Language

* R

## Bioinformatics Tools

* Seurat
* Space Ranger
* SPACEXR

## R Packages

* Seurat
* ggplot2
* dplyr
* clusterProfiler
* org.Hs.eg.db

---

# Project Structure

Spatial_Breast_Cancer/

├── scripts/

Analysis scripts for spatial transcriptomics workflow. The scripts listed
under Reproducibility form the canonical pipeline; earlier/alternative variants
were removed to keep the workflow unambiguous and remain retrievable from
earlier commits in the git history.

├── figures/

Generated visualization outputs

├── README.md

Project documentation

└── .gitignore

Git ignore configuration

---

# Data Availability

Raw Visium sequencing data, Space Ranger output files, and large Seurat objects are excluded from this repository because of size limitations.

This repository contains:

* Analysis scripts
* Visualization outputs
* Documentation
* Reproducible workflow

---

# Reproducibility

Run the analysis scripts sequentially:

1. 01_load_visium.R
2. 02_qc_filtering.R
3. 03_sctransform.R
4. 04_pca_visium.R
5. 05_clustering_visium.R
6. 06_umap_spatial_visium.R
7. 07_marker_genes_visium.R
8. 08b_annotation_validation.R
9. 09_cluster_annotation.R
10. 10_differential_expression_visium.R
11. 11_GO_enrichment_visium.R
12. 12_KEGG_enrichment_visium.R
13. 13_spatial_variable_genes.R
14. 14_plot_spatial_variable_genes.R

---

# Repository Notes

Large intermediate files are excluded:

* Seurat `.rds` objects
* Raw sequencing data
* Space Ranger outputs
* External package repositories

These files are ignored using `.gitignore`.

---

# License

Academic and research use only.


## Plain-English Explanation
See [EXPLAIN_PLAIN_ENGLISH.md](EXPLAIN_PLAIN_ENGLISH.md) for a non-technical summary of what was done, steps, and verification.
