library(Seurat)
library(ggplot2)
library(patchwork)

cat("=========================================\n")
cat(" Step 2 : Quality Control Filtering\n")
cat("=========================================\n\n")

# Create output directories
dir.create("results/QC", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Load Seurat object
breast <- readRDS("results/01_loaded_visium.rds")

cat("Object loaded successfully.\n\n")

cat("Initial dimensions:\n")
print(dim(breast))

initial_spots <- ncol(breast)

############################################################
# QC Violin Plot BEFORE filtering
############################################################

vln_before <- VlnPlot(
  breast,
  features = c(
    "nCount_Spatial",
    "nFeature_Spatial",
    "percent.mt"
  ),
  pt.size = 0.1,
  ncol = 3
)

ggsave(
  filename = "results/QC/QC_before_filtering.png",
  plot = vln_before,
  width = 12,
  height = 4,
  dpi = 300
)

############################################################
# Scatter plots
############################################################

p1 <- FeatureScatter(
  breast,
  feature1 = "nCount_Spatial",
  feature2 = "nFeature_Spatial"
)

p2 <- FeatureScatter(
  breast,
  feature1 = "nCount_Spatial",
  feature2 = "percent.mt"
)

scatter <- p1 + p2

ggsave(
  filename = "results/QC/QC_scatter.png",
  plot = scatter,
  width = 10,
  height = 5,
  dpi = 300
)

############################################################
# Filtering
############################################################

breast <- subset(
  breast,
  subset =
    nFeature_Spatial > 500 &
    nFeature_Spatial < 8000 &
    percent.mt < 20
)

filtered_spots <- ncol(breast)

cat("\nFiltering completed.\n\n")

cat("Remaining spots:", filtered_spots, "\n")
cat("Removed spots:", initial_spots - filtered_spots, "\n\n")

############################################################
# QC AFTER filtering
############################################################

vln_after <- VlnPlot(
  breast,
  features = c(
    "nCount_Spatial",
    "nFeature_Spatial",
    "percent.mt"
  ),
  pt.size = 0.1,
  ncol = 3
)

ggsave(
  filename = "results/QC/QC_after_filtering.png",
  plot = vln_after,
  width = 12,
  height = 4,
  dpi = 300
)

############################################################
# Save filtered object
############################################################

saveRDS(
  breast,
  file = "results/02_filtered_visium.rds"
)

############################################################
# Save QC summary
############################################################

qc_summary <- data.frame(
  Metric = c(
    "Initial Spots",
    "Filtered Spots",
    "Removed Spots"
  ),
  Value = c(
    initial_spots,
    filtered_spots,
    initial_spots - filtered_spots
  )
)

write.csv(
  qc_summary,
  "results/tables/QC_summary.csv",
  row.names = FALSE
)

cat("Filtered object saved.\n")
cat("results/02_filtered_visium.rds\n\n")

cat("=========================================\n")
cat(" Step 2 completed successfully\n")
cat("=========================================\n")
