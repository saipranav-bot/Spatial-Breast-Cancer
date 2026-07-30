# ==========================================
# Step 10: Differential Expression Analysis
# ==========================================

library(Seurat)
library(dplyr)

cat("====================================\n")
cat("Starting differential expression\n")
cat("====================================\n")


# Load annotated object

visium <- readRDS(
  "results/09_annotated_visium.rds"
)


# Check annotations

print(
  table(visium$cell_annotation)
)


# ------------------------------------------
# Find markers for annotated regions
# ------------------------------------------

DE_markers <- FindAllMarkers(
  visium,
  assay = "SCT",
  group.by = "cell_annotation",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


cat("Differential expression completed\n")


# Save complete DE table

write.csv(
  DE_markers,
  "results/10_DE_all_regions.csv",
  row.names = FALSE
)


# ------------------------------------------
# Top 10 genes per region
# ------------------------------------------

top10_DE <- DE_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )


write.csv(
  top10_DE,
  "results/10_DE_top10_regions.csv",
  row.names = FALSE
)


# ------------------------------------------
# Heatmap of top DE genes
# ------------------------------------------

pdf(
"figures/10_DE_heatmap.pdf",
width=10,
height=12
)

DoHeatmap(
visium,
features = unique(top10_DE$gene),
size = 2.5
) +
NoLegend()

dev.off()


cat("====================================\n")
cat("DE ANALYSIS COMPLETED\n")
cat("Saved:\n")
cat("results/10_DE_all_regions.csv\n")
cat("results/10_DE_top10_regions.csv\n")
cat("====================================\n")

