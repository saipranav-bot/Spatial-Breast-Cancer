library(Seurat)
library(ggplot2)
library(patchwork)

cat("=========================================\n")
cat(" Step 3 : SCTransform Normalization\n")
cat("=========================================\n\n")

dir.create("results/Normalization",
           recursive = TRUE,
           showWarnings = FALSE)

breast <- readRDS("results/02_filtered_visium.rds")

cat("Input object:\n")
print(breast)

############################################################
# SCTransform
############################################################

breast <- SCTransform(
    breast,
    assay = "Spatial",
    new.assay.name = "SCT",
    verbose = TRUE
)

############################################################
# Variable Feature Plot
############################################################

p <- LabelPoints(
    plot = VariableFeaturePlot(
        breast,
        assay = "SCT"
    ),
    points = head(VariableFeatures(breast),10),
    repel = TRUE
) +
labs(
    title = "Highly Variable Features",
    subtitle = "SCTransform Normalization",
    x = "Average Expression",
    y = "Standardized Variance"
) +
theme_classic(base_size = 14) +
theme(
    plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 18
    ),
    plot.subtitle = element_text(
        hjust = 0.5,
        size = 12
    ),
    axis.title = element_text(
        face = "bold",
        size = 14
    ),
    axis.text = element_text(
        colour = "black",
        size = 11
    ),
    panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.8
    ),
    legend.position = "none"
)

ggsave(
    filename = "results/Normalization/Variable_Features.png",
    plot = p,
    width = 9,
    height = 7,
    dpi = 300,
    bg = "white"
)

############################################################
# Save object
############################################################

saveRDS(
    breast,
    "results/03_sctransformed_visium.rds"
)

cat("\nSaved object:\n")
cat("results/03_sctransformed_visium.rds\n")

cat("\nSaved figure:\n")
cat("results/Normalization/Variable_Features.png\n")

cat("\n=========================================\n")
cat(" Step 3 completed successfully\n")
cat("=========================================\n")
