
library(dplyr)

cat("====================================\n")
cat("Generating cluster marker summary\n")
cat("====================================\n")


# Load top marker file

markers <- read.csv(
  "results/07_top5_visualization_markers.csv"
)


# Sort by cluster and fold change

markers <- markers %>%
  arrange(
    as.numeric(cluster),
    desc(avg_log2FC)
  )


# Print markers cluster-wise

cluster_markers <- split(
  markers$gene,
  markers$cluster
)


cat("\nTop marker genes per cluster:\n\n")


for (cl in names(cluster_markers)) {

  cat(
    "Cluster ",
    cl,
    ": ",
    paste(
      cluster_markers[[cl]],
      collapse = ", "
    ),
    "\n\n"
  )

}


# Save summary table

summary_table <- markers %>%
  group_by(cluster) %>%
  summarise(
    top_genes = paste(
      gene,
      collapse = ", "
    )
  )


write.csv(
  summary_table,
  "results/08_cluster_marker_summary.csv",
  row.names = FALSE
)


cat("====================================\n")
cat("Cluster summary completed\n")
cat("Saved: results/08_cluster_marker_summary.csv\n")
cat("====================================\n")

