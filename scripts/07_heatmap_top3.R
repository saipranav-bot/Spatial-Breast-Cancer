library(Seurat)
library(dplyr)
library(ggplot2)

cat("====================================\n")
cat("Generating Top 3 Marker Heatmap\n")
cat("====================================\n")

# ------------------------------------------
# Load object
# ------------------------------------------

visium <- readRDS("results/09_annotated_visium.rds")
DefaultAssay(visium) <- "SCT"

cat("Loaded Seurat object\n")

# ------------------------------------------
# Load marker table
# ------------------------------------------

markers <- read.csv(
  "results/07_all_cluster_markers.csv"
)

cat("Loaded", nrow(markers), "marker genes\n")

# ------------------------------------------
# Select top 3 genes per cluster
# ------------------------------------------

top3 <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 3,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  top3,
  "results/07_top3_cluster_markers.csv",
  row.names = FALSE
)

heatmap_genes <- unique(top3$gene)

cat(
  "Genes selected:",
  length(heatmap_genes),
  "\n"
)

# ------------------------------------------
# Scale only selected genes
# ------------------------------------------

visium <- ScaleData(
  visium,
  features = heatmap_genes,
  verbose = FALSE
)

cat("Scaling completed\n")

# ------------------------------------------
# Heatmap
# ------------------------------------------

pdf(
  "figures/07_top3_marker_heatmap.pdf",
  width = 16,
  height = 10
)

p <- DoHeatmap(
  visium,
  features = heatmap_genes,
  size = 4,
  raster = TRUE,
  lines.width = 0.5
) +
ggtitle("Top 3 Marker Genes per Cluster") +
theme(
  plot.title = element_text(
    size = 18,
    face = "bold",
    hjust = 0.5
  ),
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_text(
    size = 10,
    face = "bold"
  ),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 10)
)

print(p)

dev.off()

cat("Heatmap saved successfully\n")
cat("====================================\n")
cat("Step completed successfully\n")
cat("====================================\n")
