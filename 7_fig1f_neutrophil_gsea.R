library(readr)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

# ============================================================
# Figure: GSEA on DE results from TREM1+ neutrophils
# Input DE table: TET2 vs Control
# ============================================================

# ----------------------------
# paths
# ----------------------------
de_file  <- "/home/rstudio/trem1/neutrophil_v2_012026/de_cellchat/trem1pos_neut/de_output/TET2_vs_Control.tsv"
fig_dir  <- "/home/rstudio/trem1/neutrophil_v2_012026/figures"
data_dir <- "/home/rstudio/trem1/neutrophil_v2_012026/source_data"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# load DE results
# ----------------------------
df_de <- read_tsv(de_file, show_col_types = FALSE) %>%
  filter(!is.na(gene), !is.na(log2FoldChange)) %>%
  mutate(
    pvalue = ifelse(pvalue == 0, .Machine$double.xmin, pvalue),
    padj   = ifelse(padj == 0, .Machine$double.xmin, padj)
  )

write_tsv(df_de, file.path(data_dir, "panelF_de_input_source_data.tsv"))

# ----------------------------
# map SYMBOL -> ENTREZ
# ----------------------------
df_ids <- bitr(
  df_de$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
) %>%
  distinct(SYMBOL, .keep_all = TRUE)

df_kegg <- df_de %>%
  inner_join(df_ids, by = c("gene" = "SYMBOL")) %>%
  distinct(gene, .keep_all = TRUE)

# ranked vector for GSEA
ranked_gene_list <- df_kegg$log2FoldChange
names(ranked_gene_list) <- df_kegg$ENTREZID
ranked_gene_list <- sort(ranked_gene_list, decreasing = TRUE)

# ----------------------------
# run GSEA
# ----------------------------
gsea_kegg <- gseKEGG(
  geneList      = ranked_gene_list,
  organism      = "hsa",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  keyType       = "ncbi-geneid"
)

saveRDS(gsea_kegg, file.path(data_dir, "panelF_gsea_kegg_result.rds"))
write_tsv(as_tibble(gsea_kegg@result), file.path(data_dir, "panelF_gsea_kegg_table.tsv"))

# ----------------------------
# select pathway of interest
# ----------------------------
pathway_name <- "Neutrophil extracellular trap formation"

row_pathway <- gsea_kegg@result %>%
  as_tibble() %>%
  filter(Description == pathway_name)

if (nrow(row_pathway) == 0) {
  stop(paste("Pathway not found:", pathway_name))
}

pathway_id <- row_pathway$ID[1]
pathway_nes <- round(row_pathway$NES[1], 2)
pathway_p <- signif(row_pathway$pvalue[1], 3)

# ----------------------------
# plot
# ----------------------------
pdf(file.path(fig_dir, "panelF_gsea_netosis.pdf"), width = 6.5, height = 4.5)
gseaplot2(
  gsea_kegg,
  geneSetID = pathway_id,
  title = paste0(
    "GSEA: ", pathway_name, "\n",
    "NES = ", pathway_nes, ", p = ", pathway_p
  )
)
dev.off()

png(file.path(fig_dir, "panelF_gsea_netosis.png"), width = 6.5, height = 4.5, units = "in", res = 300)
gseaplot2(
  gsea_kegg,
  geneSetID = pathway_id,
  title = paste0(
    "GSEA: ", pathway_name, "\n",
    "NES = ", pathway_nes, ", p = ", pathway_p
  )
)
dev.off()