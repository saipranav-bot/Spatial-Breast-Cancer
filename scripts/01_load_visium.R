library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)

cat("=========================================\n")
cat(" Step 1 : Loading Visium dataset\n")
cat("=========================================\n\n")

data_dir <- "data/raw/outs"

cat("Checking dataset...\n\n")

stopifnot(dir.exists(data_dir))
stopifnot(file.exists(file.path(data_dir,"filtered_feature_bc_matrix.h5")))
stopifnot(dir.exists(file.path(data_dir,"spatial")))

cat("Dataset found.\n\n")

breast <- Load10X_Spatial(
    data.dir = data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "BreastCancer"
)

cat("Dataset loaded successfully.\n\n")

print(breast)

cat("\nDimensions\n")
print(dim(breast))

cat("\nCells (spots): ", ncol(breast), "\n")
cat("Genes: ", nrow(breast), "\n\n")

breast[["percent.mt"]] <- PercentageFeatureSet(
    breast,
    pattern="^MT-"
)

cat("Summary of QC metrics\n\n")

print(summary(breast$nCount_Spatial))
print(summary(breast$nFeature_Spatial))
print(summary(breast$percent.mt))

dir.create("results",showWarnings=FALSE)
dir.create("results/QC",showWarnings=FALSE,recursive=TRUE)

saveRDS(
    breast,
    file="results/01_loaded_visium.rds"
)

cat("\nSaved object:\n")
cat("results/01_loaded_visium.rds\n")

cat("\n=========================================\n")
cat(" Step 1 completed successfully\n")
cat("=========================================\n")
