# ==========================================
# Step 13: Spatially Variable Genes
# Checkpoint Version
# ==========================================

library(Seurat)
library(ggplot2)


cat("====================================\n")
cat("Starting SVG analysis\n")
cat("====================================\n")


# Load previous checkpoint

visium <- readRDS(
  "results/09_annotated_visium.rds"
)

cat("Loaded annotated object\n")



# Save checkpoint before calculation

saveRDS(
  visium,
  "results/checkpoint_13_before_svg.rds"
)

cat("Checkpoint saved: before SVG\n")



# Select genes

svg_features <- VariableFeatures(visium)

cat(
  "Number of genes for SVG:",
  length(svg_features),
  "\n"
)



# Run spatial variable gene detection

visium <- FindSpatiallyVariableFeatures(
  visium,
  assay = "SCT",
  selection.method = "markvariogram",
  features = svg_features,
  r.metric = 5,
  verbose = TRUE
)


cat("SVG calculation completed\n")



# Save immediately after calculation

saveRDS(
  visium,
  "results/checkpoint_13_after_svg.rds"
)

cat("Checkpoint saved: after SVG\n")



# Extract SVG genes

svg_genes <- SpatiallyVariableFeatures(
  visium,
  method = "markvariogram"
)


write.csv(
  data.frame(
    gene = svg_genes
  ),
  "results/13_spatial_variable_genes.csv",
  row.names = FALSE
)


cat(
  "Number of spatial genes:",
  length(svg_genes),
  "\n"
)



# Final save

saveRDS(
  visium,
  "results/13_spatial_variable_genes_visium.rds"
)


cat("====================================\n")
cat("Step 13 completed successfully\n")
cat("====================================\n")

