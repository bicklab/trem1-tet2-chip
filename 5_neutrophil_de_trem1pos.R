suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scuttle)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
})

metacell_dir <- "/home/rstudio/trem1/neutrophil_v2_012026/de_cellchat/trem1pos_neut/metacells_output"
seurat_rds   <- "/home/rstudio/trem1/neutrophil_v2_012026/rds_object/3_subsetted_trem1pos_neut.rds"
out_file     <- "/home/rstudio/trem1/neutrophil_v2_012026/de_cellchat/trem1pos_neut/de_output/TET2_vs_Control.tsv"

seu <- readRDS(seurat_rds)

control_md <- read_csv(
  file.path(metacell_dir, "control_metacells.csv"),
  show_col_types = FALSE
) %>%
  rename(cell_barcode = 1) %>%
  mutate(condition = "Control")

tet2_md <- read_csv(
  file.path(metacell_dir, "tet2_metacells.csv"),
  show_col_types = FALSE
) %>%
  rename(cell_barcode = 1) %>%
  mutate(condition = "TET2")

md <- bind_rows(control_md, tet2_md) %>%
  filter(metacell >= 0)

counts <- GetAssayData(seu, assay = "RNA", slot = "counts")

md <- md %>%
  filter(cell_barcode %in% colnames(counts)) %>%
  mutate(
    condition = factor(condition, levels = c("Control", "TET2")),
    metacell_id = paste(condition, metacell, sep = "__")
  )

counts_sub <- counts[, md$cell_barcode, drop = FALSE]

sce <- SingleCellExperiment(
  assays = list(counts = counts_sub),
  colData = md
)

agg <- aggregateAcrossCells(
  sce,
  ids = DataFrame(
    condition = colData(sce)$condition,
    metacell_id = colData(sce)$metacell_id
  )
)

coldata <- as.data.frame(colData(agg))
coldata$condition <- factor(coldata$condition, levels = c("Control", "TET2"))

dds <- DESeqDataSetFromMatrix(
  countData = round(assay(agg, "counts")),
  colData = coldata,
  design = ~ condition
)

keep <- rowSums(counts(dds) >= 6) >= ceiling(0.6 * ncol(dds))
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "TET2", "Control"))

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_tsv(res_df, out_file)