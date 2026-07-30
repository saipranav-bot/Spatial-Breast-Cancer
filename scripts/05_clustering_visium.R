# ==========================================
# Step 05: Find Neighbors and Clustering
# ==========================================

library(Seurat)
library(ggplot2)

cat("====================================\n")
cat("Starting clustering analysis\n")
cat("====================================\n")


# Load PCA object
visium <- readRDS(
  "results/04_pca_visium.rds"
)

cat("Loaded PCA object\n")


# Number of PCs selected from elbow plot
pcs.use <- 1:20


# Find nearest neighbors
visium <- FindNeighbors(
  visium,
  reduction = "pca",
  dims = pcs.use,
  verbose = TRUE
)

cat("Neighbor graph completed\n")


# Cluster spots
visium <- FindClusters(
  visium,
  resolution = 0.5,
  verbose = TRUE
)

cat("Clustering completed\n")


# Cluster plot
pdf(
  "figures/05_cluster_pca.pdf",
  width = 8,
  height = 6
)

DimPlot(
  visium,
  reduction = "pca",
  group.by = "seurat_clusters",
  label = TRUE
)

dev.off()


# Save clustered object
saveRDS(
  visium,
  "results/05_clustered_visium.rds"
)


cat("====================================\n")
cat("DONE\n")
cat("Saved: results/05_clustered_visium.rds\n")
cat("====================================\n")
