# ==========================================
# Step 06: UMAP and Spatial Cluster Plot
# ==========================================

library(Seurat)
library(ggplot2)

cat("====================================\n")
cat("Starting UMAP analysis\n")
cat("====================================\n")


# Load clustered object
visium <- readRDS(
  "results/05_clustered_visium.rds"
)

cat("Loaded clustered Visium object\n")


# -----------------------------
# Run UMAP
# -----------------------------

visium <- RunUMAP(
  visium,
  reduction = "pca",
  dims = 1:20,
  verbose = TRUE
)

cat("UMAP completed\n")


# -----------------------------
# UMAP cluster plot
# -----------------------------

pdf(
  "figures/06_umap_clusters.pdf",
  width = 8,
  height = 6
)

DimPlot(
  visium,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

dev.off()

cat("UMAP plot saved\n")


# -----------------------------
# Spatial cluster plot
# -----------------------------

pdf(
  "figures/06_spatial_clusters.pdf",
  width = 8,
  height = 8
)

SpatialDimPlot(
  visium,
  group.by = "seurat_clusters",
  label = TRUE
)

dev.off()

cat("Spatial cluster plot saved\n")


# -----------------------------
# Save UMAP object
# -----------------------------

saveRDS(
  visium,
  "results/06_umap_visium.rds"
)


cat("====================================\n")
cat("UMAP + SPATIAL CLUSTERING COMPLETED\n")
cat("Saved: results/06_umap_visium.rds\n")
cat("====================================\n")

