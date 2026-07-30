library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(jsonlite)


data_dir <- "data/raw/V1_Human_Breast_Cancer"

matrix_dir <- file.path(
  data_dir,
  "filtered_feature_bc_matrix"
)

spatial_dir <- file.path(
  data_dir,
  "spatial"
)


# Read expression matrix

counts <- Read10X(
  data.dir = matrix_dir
)


# Create Seurat object

breast <- CreateSeuratObject(
  counts = counts,
  assay = "Spatial",
  project = "Human_Breast_Cancer"
)

print(breast)


# Add mitochondrial percentage

breast[["percent.mt"]] <- PercentageFeatureSet(
  breast,
  pattern = "^MT-"
)


# QC summary

print(head(breast[[]]))


# Filter spots

breast <- subset(
  breast,
  subset =
    nFeature_Spatial > 500 &
    nFeature_Spatial < 8000 &
    percent.mt < 20
)

print(breast)


# QC plot

p <- VlnPlot(
  breast,
  features = c(
    "nCount_Spatial",
    "nFeature_Spatial",
    "percent.mt"
  ),
  pt.size = 0.1,
  ncol = 3
)

ggsave(
  "results/figures/QC_after_filtering.png",
  p,
  width = 12,
  height = 4
)


# SCTransform normalization

breast <- SCTransform(
  breast,
  assay = "Spatial",
  verbose = FALSE
)


# PCA

breast <- RunPCA(
  breast,
  assay = "SCT",
  verbose = FALSE
)


# Save object

saveRDS(
  breast,
  "results/breast_visium_sct.rds"
)


cat("QC normalization PCA completed successfully\n")
