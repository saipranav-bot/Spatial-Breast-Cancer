library(Seurat)
library(dplyr)
library(pheatmap)

cat("====================================\n")
cat("Cluster Average Heatmap\n")
cat("====================================\n")

# Load object
visium <- readRDS("results/09_annotated_visium.rds")
DefaultAssay(visium) <- "SCT"

# Read markers
markers <- read.csv("results/07_all_cluster_markers.csv")

# Top 3 markers per cluster
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

genes <- unique(top3$gene)

cat("Genes selected:", length(genes), "\n")

# Average expression for each cluster
avg <- AverageExpression(
  visium,
  assays = "SCT",
  features = genes,
  group.by = "seurat_clusters",
  return.seurat = FALSE
)

mat <- avg$SCT

# Z-score each gene
mat <- t(scale(t(mat)))

pdf(
  "figures/07_cluster_average_heatmap.pdf",
  width = 8,
  height = 10
)

pheatmap(
  mat,
  color = colorRampPalette(
    c("navy","white","firebrick3")
  )(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 12,
  angle_col = 0,
  border_color = "grey80",
  main = "Top 3 Marker Genes per Cluster"
)

dev.off()

cat("Heatmap saved successfully\n")
cat("====================================\n")
cat("Completed\n")
cat("====================================\n")
