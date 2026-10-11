library(Seurat)
library(ggplot2)
library(patchwork)

# Load the QC-filtered object
breast <- readRDS("results/breast_visium_sct.rds")

print(breast)

# -----------------------
# PCA
# -----------------------

breast <- RunPCA(
  breast,
  assay = "SCT",
  verbose = FALSE
)

# Elbow Plot
png(
  "results/figures/PCA_Elbow.png",
  width = 900,
  height = 700
)

ElbowPlot(breast, ndims = 50)

dev.off()

# -----------------------
# Neighbor graph
# -----------------------

breast <- FindNeighbors(
  breast,
  reduction = "pca",
  dims = 1:30
)

# -----------------------
# Clustering
# -----------------------

breast <- FindClusters(
  breast,
  resolution = 0.5
)

# -----------------------
# UMAP
# -----------------------

breast <- RunUMAP(
  breast,
  reduction = "pca",
  dims = 1:30
)

# Save UMAP

p1 <- DimPlot(
  breast,
  reduction = "umap",
  label = TRUE
)

ggsave(
  "results/figures/UMAP_clusters.png",
  p1,
  width = 8,
  height = 6
)

# Save PCA

p2 <- DimPlot(
  breast,
  reduction = "pca"
)

ggsave(
  "results/figures/PCA.png",
  p2,
  width = 8,
  height = 6
)

saveRDS(
  breast,
  "results/breast_visium_clustered.rds"
)

cat("\n====================================\n")
cat("Dimensionality reduction completed!\n")
cat("====================================\n")
cat("Number of clusters:\n")
print(table(Idents(breast)))
