library(Seurat)
library(ggplot2)

cat("=========================================\n")
cat(" Step 3 : SCTransform Normalization\n")
cat("=========================================\n\n")

dir.create("results/Normalization", recursive = TRUE, showWarnings = FALSE)

#--------------------------------------------------
# Load filtered object
#--------------------------------------------------

breast <- readRDS("results/02_filtered_visium.rds")

cat("Dataset loaded.\n\n")

print(breast)

cat("\nRunning SCTransform...\n\n")

#--------------------------------------------------
# SCTransform
#--------------------------------------------------

breast <- SCTransform(
    breast,
    assay = "Spatial",
    new.assay.name = "SCT",
    variable.features.n = 3000,
    verbose = TRUE
)

cat("\nSCTransform completed.\n\n")

#--------------------------------------------------
# Summary
#--------------------------------------------------

cat("Active assay:\n")
print(DefaultAssay(breast))

cat("\nVariable features:\n")
print(length(VariableFeatures(breast)))

cat("\nTop 20 variable genes:\n")
print(head(VariableFeatures(breast),20))

#--------------------------------------------------
# Variable feature plot
#--------------------------------------------------

p <- VariableFeaturePlot(breast)

ggsave(
    filename="results/Normalization/Variable_Features.png",
    plot=p,
    width=8,
    height=6,
    dpi=300
)

#--------------------------------------------------
# Save object
#--------------------------------------------------

saveRDS(
    breast,
    file="results/03_sctransform.rds"
)

cat("\nSaved object:\n")
cat("results/03_sctransform.rds\n")

cat("\n=========================================\n")
cat(" Step 3 completed successfully\n")
cat("=========================================\n")
