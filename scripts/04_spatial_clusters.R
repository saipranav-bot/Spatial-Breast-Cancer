library(Seurat)
library(ggplot2)

# Load clustered object
breast <- readRDS("results/breast_visium_clustered.rds")

print(breast)

# Plot clusters on tissue
p1 <- SpatialDimPlot(
  breast,
  label = TRUE,
  label.size = 4,
  pt.size.factor = 1.6
)

ggsave(
  filename = "results/figures/Spatial_clusters.png",
  plot = p1,
  width = 8,
  height = 8,
  dpi = 300
)

cat("\n=====================================\n")
cat("Spatial cluster plot saved!\n")
cat("=====================================\n")
