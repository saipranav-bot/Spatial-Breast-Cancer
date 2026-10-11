# ==========================================
# Publication-quality Marker Gene Plots
# ==========================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("====================================\n")
cat("Publication Marker Plot Pipeline\n")
cat("====================================\n")

# ------------------------------------------
# Load Seurat object
# ------------------------------------------

visium <- readRDS("results/09_annotated_visium.rds")
DefaultAssay(visium) <- "SCT"

cat("Loaded Seurat object\n")

# ------------------------------------------
# Load marker table
# ------------------------------------------

markers <- read.csv("results/07_all_cluster_markers.csv")

cat("Loaded", nrow(markers), "marker genes\n")

# ------------------------------------------
# Choose ONE UNIQUE marker per cluster
# ------------------------------------------

candidate <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  ungroup()

selected <- data.frame()

used.genes <- character()

for(cl in sort(unique(candidate$cluster))){

    tmp <- candidate %>%
      filter(cluster == cl)

    gene <- tmp$gene[!tmp$gene %in% used.genes][1]

    selected <- rbind(
      selected,
      tmp[tmp$gene == gene, ]
    )

    used.genes <- c(used.genes, gene)

}

genes <- selected$gene

cat("Representative genes:\n")
print(genes)

write.csv(
  selected,
  "results/07_top1_cluster_markers.csv",
  row.names = FALSE
)

# ------------------------------------------
# Scale selected genes
# ------------------------------------------

visium <- ScaleData(
  visium,
  features = genes,
  verbose = FALSE
)

cat("Scaling completed\n")

# ------------------------------------------
# Heatmap
# ------------------------------------------

pdf(
  "figures/07_publication_heatmap.pdf",
  width = 10,
  height = 8
)

p1 <- DoHeatmap(
  visium,
  features = genes,
  group.by = "seurat_clusters",
  raster = TRUE
) +
theme_classic(base_size = 14) +
theme(
  axis.text.y = element_text(size = 11, face = "bold"),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 16
  )
) +
ggtitle("Representative Marker Gene per Cluster")

print(p1)

dev.off()

cat("Heatmap saved\n")

# ------------------------------------------
# DotPlot
# ------------------------------------------

pdf(
  "figures/07_publication_dotplot.pdf",
  width = 9,
  height = 6
)

p2 <- DotPlot(
  visium,
  features = genes
) +
RotatedAxis() +
theme_classic(base_size = 14) +
theme(
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    size = 10,
    face = "bold"
  ),
  axis.text.y = element_text(
    size = 11,
    face = "bold"
  ),
  plot.title = element_text(
    hjust = 0.5,
    size = 16,
    face = "bold"
  )
) +
ggtitle("Representative Marker Genes")

print(p2)

dev.off()

cat("DotPlot saved\n")

# ------------------------------------------
# Spatial Feature Plots
# ------------------------------------------

pdf(
  "figures/07_top_marker_spatial.pdf",
  width = 16,
  height = 12
)

plots <- lapply(
  genes,
  function(g){

    SpatialFeaturePlot(
      visium,
      features = g
    ) +
    ggtitle(g) +
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
    )

  }
)

print(
  wrap_plots(
    plots,
    ncol = 3
  )
)

dev.off()

cat("Spatial plots saved\n")

cat("====================================\n")
cat("Pipeline completed successfully\n")
cat("====================================\n")
