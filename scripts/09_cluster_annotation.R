# ==========================================
# Step 09: Spatial Cluster Annotation
# ==========================================

library(Seurat)

cat("====================================\n")
cat("Adding cluster annotations\n")
cat("====================================\n")


# Load object

visium <- readRDS(
  "results/07_marker_visium.rds"
)


# Current clusters

print(
  table(visium$seurat_clusters)
)


# Create annotation mapping

cluster_labels <- c(
  
  "0" = "TLS_Lymphoid_Immune",
  "1" = "Ciliated_Epithelial",
  "2" = "Hypoxic_TAM_Stromal",
  "3" = "Basal_Tumor_Epithelial",
  "4" = "Adipocyte_Stromal",
  "5" = "Luminal_Epithelial",
  "6" = "Neural_like_Stromal",
  "7" = "Secretory_Epithelial",
  "8" = "Immune_Associated",
  "9" = "Differentiated_Epithelial",
  "10" = "Neural_like_Region",
  "11" = "Secretory_Epithelial_2"
  
)


# Add annotation column

visium$cell_annotation <- plyr::mapvalues(
  visium$seurat_clusters,
  from = names(cluster_labels),
  to = cluster_labels
)


# Check

table(visium$cell_annotation)


# Save annotated object

saveRDS(
  visium,
  "results/09_annotated_visium.rds"
)


# UMAP annotation plot

pdf(
"figures/09_annotated_umap.pdf",
width=10,
height=8
)

DimPlot(
visium,
reduction="umap",
group.by="cell_annotation",
label=TRUE
)

dev.off()


# Spatial annotation plot

pdf(
"figures/09_annotated_spatial.pdf",
width=10,
height=10
)

SpatialDimPlot(
visium,
group.by="cell_annotation",
label=TRUE
)

dev.off()


cat("====================================\n")
cat("Annotation completed\n")
cat("Saved: results/09_annotated_visium.rds\n")
cat("====================================\n")

