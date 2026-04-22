# ============================================================
# Preprocessing and QC for whole-blood neutrophil analysis
# 1) remove ambient RNA with decontX
# 2) filter low-quality cells
# 3) detect doublets with DoubletFinder
# 4) re-normalize singlets for clustering/annotation
# 5) annotate cell types with ScType
# 6) save processed Seurat object
# ============================================================

.libPaths("/home/rstudio/4.4-3.20")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SingleCellExperiment)
  library(celda)
  library(DoubletFinder)
})

options(future.globals.maxSize = 4 * 1024^3)

input_dir  <- "/home/rstudio/trem1/neutrophil_v2_012026/input/3p_10x_wholeblood_control_matrix/"
output_rds <- "/home/rstudio/trem1/neutrophil_v2_012026/rds_object/1_annotated_3p_wb_control_rnacount_800.rds"

remove_ambient <- function(input_dir) {
  counts <- Read10X(data.dir = input_dir)
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- decontX(sce)
  CreateSeuratObject(counts = round(decontXcounts(sce)))
}

quality_control <- function(seu) {
  seu <- subset(seu, subset = nFeature_RNA >= 200 & nCount_RNA >= 800)
  
  seu <- PercentageFeatureSet(
    seu,
    pattern = "^MT-",
    col.name = "percent.mt",
    assay = "RNA"
  )
  
  seu <- subset(seu, subset = percent.mt < 25)
  
  rm_genes <- rownames(seu)[grepl("^MT-|^MRP", rownames(seu))]
  subset(seu, features = setdiff(rownames(seu), rm_genes))
}

run_sct_clustering <- function(seu, dims = 1:10) {
  seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = dims, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = dims, verbose = FALSE)
  seu <- FindClusters(seu, verbose = FALSE)
  seu
}

find_pK_noplot <- function(sweep_stats) {
  if ("AUC" %in% colnames(sweep_stats)) {
    stop("AUC mode not supported in this version.")
  }
  
  pks <- unique(sweep_stats$pK)
  out <- data.frame(
    pK = pks,
    MeanBC = NA_real_,
    VarBC = NA_real_,
    BCmetric = NA_real_
  )
  
  for (i in seq_along(pks)) {
    idx <- which(sweep_stats$pK == pks[i])
    mean_bc <- mean(sweep_stats[idx, "BCreal"])
    sd_bc <- sd(sweep_stats[idx, "BCreal"])
    out$MeanBC[i] <- mean_bc
    out$VarBC[i] <- sd_bc^2
    out$BCmetric[i] <- mean_bc / (sd_bc^2)
  }
  
  out
}

run_doublet_finder <- function(seu, pcs_use = 1:10) {
  sweep_res <- paramSweep(seu, PCs = pcs_use, sct = TRUE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find_pK_noplot(sweep_stats)
  
  best_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  best_pK <- as.numeric(as.character(best_pK))
  
  homotypic.prop <- modelHomotypic(seu$seurat_clusters)
  doublet_rate <- (ncol(seu) / 0.57) * 4.6e-6
  nExp <- round(doublet_rate * ncol(seu))
  nExp_adj <- round(nExp * (1 - homotypic.prop))
  
  seu <- doubletFinder(
    seu,
    PCs = pcs_use,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp_adj,
    sct = TRUE
  )
  
  df_col <- colnames(seu@meta.data)[str_detect(colnames(seu@meta.data), "DF.classifications")]
  pann_col <- colnames(seu@meta.data)[str_detect(colnames(seu@meta.data), "^pANN")]
  
  seu$DoubletFinder <- seu@meta.data[[df_col]]
  seu@meta.data[[df_col]] <- NULL
  
  seu$pANN <- seu@meta.data[[pann_col]]
  seu@meta.data[[pann_col]] <- NULL
  
  seu
}

prepare_for_annotation <- function(seu) {
  seu <- subset(seu, subset = DoubletFinder == "Singlet")
  
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
  
  seu <- SCTransform(
    seu,
    assay = "RNA",
    new.assay.name = "SCT",
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
  
  DefaultAssay(seu) <- "SCT"
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:15, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:15, verbose = FALSE)
  seu <- FindClusters(seu, verbose = FALSE)
  
  rm_genes <- rownames(seu)[grepl("^RPS|^RPL|^MT-|^MRP", rownames(seu))]
  subset(seu, features = setdiff(rownames(seu), rm_genes))
}

annotate_with_sctype <- function(seu) {
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  gs_list <- gene_sets_prepare(
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
    "Immune system"
  )
  
  es.max <- sctype_score(
    scRNAseqData = seu[["RNA"]]@scale.data,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  
  cluster_scores <- do.call(
    rbind,
    lapply(unique(seu$seurat_clusters), function(cl) {
      cl_cells <- rownames(seu@meta.data[seu$seurat_clusters == cl, , drop = FALSE])
      scores <- sort(rowSums(es.max[, cl_cells, drop = FALSE]), decreasing = TRUE)
      head(
        data.frame(
          cluster = cl,
          type = names(scores),
          scores = scores,
          ncells = length(cl_cells)
        ),
        10
      )
    })
  )
  
  top_calls <- cluster_scores %>%
    group_by(cluster) %>%
    slice_max(order_by = scores, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  top_calls$type[as.numeric(top_calls$scores) < top_calls$ncells / 10] <- "Unknown"
  
  seu$customclassif <- ""
  for (cl in unique(top_calls$cluster)) {
    seu$customclassif[seu$seurat_clusters == cl] <- as.character(top_calls$type[top_calls$cluster == cl][1])
  }
  
  seu
}

message("Step 1: ambient RNA removal")
seu <- remove_ambient(input_dir)

message("Step 2: cell-level QC")
seu <- quality_control(seu)

message("Step 3: clustering for DoubletFinder")
seu <- run_sct_clustering(seu, dims = 1:10)

message("Step 4: doublet calling")
seu <- run_doublet_finder(seu, pcs_use = 1:10)

message("Step 5: reprocessing singlets for annotation")
seu <- prepare_for_annotation(seu)

message("Step 6: cell-type annotation")
seu <- annotate_with_sctype(seu)

# manual relabeling used in downstream analysis
seu$customclassif <- as.character(seu$customclassif)
seu$customclassif[seu$customclassif == "Unknown"] <- "Neutrophils"
seu$customclassif[seu$customclassif == "Myeloid Dendritic cells"] <- "Non-classical Monocytes"

message("Saving processed object")
saveRDS(seu, output_rds)

message("Done")