# ==========================================
# Step 14: Plot Top 5 Spatially Variable Genes
# ==========================================

library(Seurat)
library(ggplot2)
library(patchwork)

cat("====================================\n")
cat("Plotting Top 5 Spatially Variable Genes\n")
cat("====================================\n")


# ------------------------------------------
# Load SVG object
# ------------------------------------------

visium <- readRDS(
  "results/13_spatial_variable_genes_visium.rds"
)

cat("SVG object loaded\n")


# ------------------------------------------
# Get SVG genes
# ------------------------------------------

svg_genes <- SpatiallyVariableFeatures(
  visium,
  method = "markvariogram"
)

cat(
  "Number of spatially variable genes:",
  length(svg_genes),
  "\n"
)


# ------------------------------------------
# Save complete SVG list
# ------------------------------------------

write.csv(
  data.frame(Gene = svg_genes),
  "results/14_spatial_variable_gene_list.csv",
  row.names = FALSE
)


# ------------------------------------------
# Select Top 5 genes
# ------------------------------------------

top5 <- head(svg_genes, 5)

write.csv(
  data.frame(Gene = top5),
  "results/14_top5_spatial_variable_genes.csv",
  row.names = FALSE
)

cat("\nTop 5 Spatially Variable Genes:\n")
print(top5)


# ------------------------------------------
# Create individual plots
# ------------------------------------------

plots <- lapply(top5, function(gene){

  SpatialFeaturePlot(
    visium,
    features = gene,
    pt.size.factor = 1.8,
    alpha = c(0.2, 1)
  ) +
    ggtitle(gene) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5
      ),
      axis.title.x = element_text(
        size = 14,
        face = "bold"
      ),
      axis.title.y = element_text(
        size = 14,
        face = "bold"
      ),
      axis.text.x = element_text(
        size = 11,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 11,
        colour = "black"
      ),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )

})


# ------------------------------------------
# Save publication-quality figure
# ------------------------------------------

pdf(
  "figures/14_top5_spatial_variable_genes.pdf",
  width = 14,
  height = 16
)

wrap_plots(
  plots,
  ncol = 1
)

dev.off()

cat("Top 5 SVG plot saved\n")
cat("====================================\n")
cat("Step 14 completed successfully\n")
cat("====================================\n")
