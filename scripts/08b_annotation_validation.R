library(Seurat)
library(ggplot2)

visium <- readRDS(
  "results/07_marker_visium.rds"
)

markers_check <- c(
"EPCAM",
"KRT8",
"KRT18",
"KRT19",
"KRT5",
"KRT14",
"TP63",
"ESR1",
"PGR",
"FOXA1",
"GATA3",
"SCGB2A2",
"CD3D",
"CD3E",
"CD79A",
"MS4A1",
"LST1",
"C1QA",
"COL1A1",
"COL3A1",
"DCN",
"LUM"
)


markers_check <- markers_check[
  markers_check %in% rownames(visium)
]


pdf(
"figures/08b_annotation_validation_dotplot.pdf",
width=12,
height=8
)

DotPlot(
visium,
features = markers_check
) +
RotatedAxis()

dev.off()


pdf(
"figures/08b_annotation_spatial_markers.pdf",
width=12,
height=12
)

SpatialFeaturePlot(
visium,
features = markers_check[1:12],
ncol=3
)

dev.off()


cat("Annotation validation plots completed\n")

