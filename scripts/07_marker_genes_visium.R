# ==========================================
# Step 07: Marker Gene Analysis
# Publication Version
# ==========================================

library(Seurat)
library(dplyr)
library(ggplot2)

cat("====================================\n")
cat("Step 07 : Marker Gene Analysis\n")
cat("====================================\n")


# ------------------------------------------
# Load annotated object
# ------------------------------------------

visium <- readRDS(
  "results/09_annotated_visium.rds"
)

DefaultAssay(visium) <- "SCT"

cat("Object loaded successfully\n")


# ------------------------------------------
# Find marker genes
# ------------------------------------------

markers <- FindAllMarkers(
  visium,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cat("Marker discovery completed\n")


write.csv(
  markers,
  "results/07_all_cluster_markers.csv",
  row.names = FALSE
)


# ------------------------------------------
# Top 10 markers per cluster
# ------------------------------------------

top10_cluster <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  )

write.csv(
  top10_cluster,
  "results/07_top10_cluster_markers.csv",
  row.names = FALSE
)


# ------------------------------------------
# Top 10 markers overall
# ------------------------------------------

top10 <- markers %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene,.keep_all=TRUE) %>%
  slice_head(n=10)

write.csv(
  top10,
  "results/07_top10_overall_markers.csv",
  row.names = FALSE
)

cat("Top 10 overall marker genes selected\n")


# ------------------------------------------
# Heatmap
# ------------------------------------------

pdf(
  "figures/07_top10_marker_heatmap.pdf",
  width=8,
  height=7
)

DoHeatmap(
  visium,
  features=top10$gene,
  size=5
)+
theme(
axis.text.x=element_text(
size=10,
angle=45,
hjust=1,
face="bold"
),
axis.text.y=element_text(
size=10,
face="bold"
),
plot.title=element_text(
size=14,
face="bold",
hjust=0.5
)
)

dev.off()


# ------------------------------------------
# DotPlot
# ------------------------------------------

pdf(
  "figures/07_top10_marker_dotplot.pdf",
  width=8,
  height=6
)

DotPlot(
  visium,
  features=top10$gene
)+
RotatedAxis()+
theme(
axis.text.x=element_text(
size=10,
angle=45,
hjust=1
),
axis.text.y=element_text(
size=10
),
plot.title=element_text(
size=14,
face="bold"
)
)

dev.off()


cat("====================================\n")
cat("Step 07 completed successfully\n")
cat("====================================\n")
