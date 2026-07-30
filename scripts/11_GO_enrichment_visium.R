# ==========================================
# Step 11: GO Enrichment (Clean Top 5 Plot)
# ==========================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)


cat("====================================\n")
cat("Starting GO enrichment analysis\n")
cat("====================================\n")


DE_markers <- read.csv(
  "results/10_DE_all_regions.csv"
)


dir.create(
  "results/GO",
  showWarnings = FALSE
)

dir.create(
  "figures/GO",
  showWarnings = FALSE
)


regions <- unique(DE_markers$cluster)


for(region in regions){

  cat("Processing:", region, "\n")


  genes <- DE_markers %>%
    filter(
      cluster == region,
      p_val_adj < 0.05,
      avg_log2FC > 0.25
    ) %>%
    pull(gene)


  gene_ids <- bitr(
    genes,
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )


  if(nrow(gene_ids) < 5){
    next
  }


  ego <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )


  ego_df <- as.data.frame(ego)


  write.csv(
    ego_df,
    paste0(
      "results/GO/",
      region,
      "_GO_BP.csv"
    ),
    row.names=FALSE
  )


  if(nrow(ego_df) > 0){

    pdf(
      paste0(
        "figures/GO/",
        region,
        "_GO_dotplot_clean.pdf"
      ),
      width=10,
      height=6
    )


    print(
      dotplot(
        ego,
        showCategory = 5,
        font.size = 10
      ) +
      ggtitle(region)
    )


    dev.off()

  }

}


cat("====================================\n")
cat("CLEAN GO PLOTS COMPLETED\n")
cat("Top 5 pathways only\n")
cat("====================================\n")

