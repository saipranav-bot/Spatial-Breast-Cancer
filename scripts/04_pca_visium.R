
# ==========================================
# Step 04: PCA Analysis - Visium Spatial Data
# ==========================================

library(Seurat)
library(ggplot2)

set.seed(123)

cat("====================================\n")
cat("Starting PCA analysis\n")
cat("====================================\n")


# Fixed seed for reproducibility
set.seed(123)

# Fixed seed for reproducibility
set.seed(123)

# Load SCT normalized object
visium <- readRDS(
  "results/03_sctransform.rds"
)

cat("Loaded SCT normalized Visium object\n")


# Run PCA
visium <- RunPCA(
  visium,
  assay = "SCT",
  npcs = 50,
  verbose = TRUE
)

cat("PCA completed successfully\n")


# Print PCA information
print(visium[["pca"]])


# Save Elbow Plot
pdf(
  "figures/04_elbow_plot.pdf",
  width = 8,
  height = 6
)

ElbowPlot(
  visium,
  ndims = 50
)

dev.off()

cat("Elbow plot saved\n")


# Save PCA plot
pdf(
  "figures/04_pca_plot.pdf",
  width = 8,
  height = 6
)

DimPlot(
  visium,
  reduction = "pca",
  dims = c(1,2)
)

dev.off()

cat("PCA plot saved\n")


# Save PCA object
saveRDS(
  visium,
  "results/04_pca_visium.rds"
)


cat("====================================\n")
cat("PCA ANALYSIS COMPLETED\n")
cat("Saved: results/04_pca_visium.rds\n")
cat("====================================\n")

