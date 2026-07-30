# ==========================================
# Step 12: KEGG Enrichment Analysis
# ==========================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)


cat("====================================\n")
cat("Starting KEGG enrichment analysis\n")
cat("====================================\n")


# Load DE markers

DE_markers <- read.csv(
  "results/10_DE_all_regions.csv"
)


# Output folders

dir.create(
  "results/KEGG",
  showWarnings = FALSE
)

dir.create(
  "figures/KEGG",
  showWarnings = FALSE
)


regions <- unique(DE_markers$cluster)


for(region in regions){

  cat("Processing:", region, "\n")


  # Select significant genes

  genes <- DE_markers %>%
    filter(
      cluster == region,
      p_val_adj < 0.05,
      avg_log2FC > 0.25
    ) %>%
    pull(gene)


  # Convert gene symbols

  gene_ids <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )


  if(nrow(gene_ids) < 5){
    cat("Skipping", region, "\n")
    next
  }


  # KEGG enrichment

  ekegg <- enrichKEGG(
    gene = gene_ids$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )


  if(is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0){
    next
  }


  # Convert IDs to symbols

  ekegg <- setReadable(
    ekegg,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )


  # Save table

  write.csv(
    as.data.frame(ekegg),
    paste0(
      "results/KEGG/",
      region,
      "_KEGG.csv"
    ),
    row.names = FALSE
  )


  # Clean top 5 plot

  pdf(
    paste0(
      "figures/KEGG/",
      region,
      "_KEGG_dotplot.pdf"
    ),
    width = 10,
    height = 6
  )


  print(
    dotplot(
      ekegg,
      showCategory = 5,
      font.size = 10
    ) +
    ggtitle(region)
  )


  dev.off()


}


cat("====================================\n")
cat("KEGG ENRICHMENT COMPLETED\n")
cat("Top 5 pathways plotted\n")
cat("====================================\n")

