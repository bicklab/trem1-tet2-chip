library(Seurat)
library(CellChat)

# ============================================================
# Figure: Mo-Mac 2 -> endothelial CellChat comparison
# TET2 vs Control in adipose tissue
# ============================================================

# ----------------------------
# 1. Load input Seurat object
# ----------------------------
# This should be the Seurat object containing the cell populations
# used for this figure, including:
# - Adipose Tissue Mo-Mac 2
# - Adipose Tissue Capillary EC
# - Adipose Tissue Venous EC
# - Adipose Tissue Arterial EC
# and a metadata column named "CHIP" with values such as "Control" and "TET2"

seu_momac2_ec_input <- readRDS(
  "/home/rstudio/trem1/human_adi/rds_obj/6_ec_adi_myeloid_trem1_pos.rds"
)

# ----------------------------
# 2. Keep only the cell types used in this figure (mo-mac2 and all endothelial cells)
# ----------------------------
celltypes_momac2_ec <- c(
  "Adipose Tissue Mo-Mac 2",
  "Adipose Tissue Capillary EC",
  "Adipose Tissue Venous EC",
  "Adipose Tissue Arterial EC"
)

seu_momac2_ec <- subset(
  x = seu_momac2_ec_input,
  subset = Cell_Type %in% celltypes_momac2_ec
)

DefaultAssay(seu_momac2_ec) <- "RNA"

# ----------------------------
# 3. Split into Control and TET2
# ----------------------------
seu_momac2_ec_ctrl <- subset(
  x = seu_momac2_ec,
  subset = CHIP == "Control"
)

seu_momac2_ec_tet2 <- subset(
  x = seu_momac2_ec,
  subset = CHIP == "TET2"
)

# ----------------------------
# 4. Build CellChat objects
# ----------------------------
build_cellchat_object <- function(seu_obj, group_by, output_rds = NULL) {
  cellchat_obj <- createCellChat(
    object = seu_obj,
    meta = seu_obj@meta.data,
    group.by = group_by
  )
  
  cellchat_obj@DB <- CellChatDB.human
  
  cellchat_obj <- subsetData(cellchat_obj)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- projectData(cellchat_obj, PPI.human)
  cellchat_obj <- computeCommunProb(
    cellchat_obj,
    type = "truncatedMean",
    trim = 0.1
  )
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 25)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(
    cellchat_obj,
    slot.name = "netP"
  )
  
  if (!is.null(output_rds)) {
    saveRDS(cellchat_obj, output_rds)
  }
  
  return(cellchat_obj)
}

cellchat_momac2_ec_ctrl <- build_cellchat_object(
  seu_obj = seu_momac2_ec_ctrl,
  group_by = "Cell_Type",
  output_rds = "/home/rstudio/trem1/human_adi/rds_obj/cellchat_momac2_ec_ctrl.rds"
)

cellchat_momac2_ec_tet2 <- build_cellchat_object(
  seu_obj = seu_momac2_ec_tet2,
  group_by = "Cell_Type",
  output_rds = "/home/rstudio/trem1/human_adi/rds_obj/cellchat_momac2_ec_tet2.rds"
)

# ----------------------------
# 5. Merge CellChat objects
# ----------------------------
cellchat_momac2_ec_merged <- mergeCellChat(
  object.list = list(
    CTRL = cellchat_momac2_ec_ctrl,
    TET2 = cellchat_momac2_ec_tet2
  ),
  add.names = c("CTRL", "TET2"),
  cell.prefix = TRUE
)

cellchat_momac2_ec_merged <- updateCellChat(cellchat_momac2_ec_merged)

# ----------------------------
# 6. Differential interaction circle plot
# ----------------------------
groups_momac2_ec <- c(
  "Adipose Tissue Mo-Mac 2",
  "Adipose Tissue Capillary EC",
  "Adipose Tissue Venous EC",
  "Adipose Tissue Arterial EC"
)

svg(
  "/home/rstudio/trem1/human_adi/cellchat_figure/cellchat_diff_circle_momac2_ec.svg",
  height = 6,
  width = 6
)

netVisual_diffInteraction(
  cellchat_momac2_ec_merged,
  weight.scale = TRUE,
  measure = "weight",
  sources.use = groups_momac2_ec,
  targets.use = groups_momac2_ec,
  remove.isolate = TRUE,
  color.use = c("#B095C6", "#8AA6D8", "#4BB"),
  color.edge = c("#72B28A", "#D88A8A")
)

dev.off()

# ----------------------------
# 7. Pathway-level comparison plot
# ----------------------------
pathways_momac2_to_ec <- c(
  "CD46",
  "OSM",
  "ITGB2",
  "CXCL",
  "SELPLG",
  "VCAM",
  "THBS",
  "RESISTIN",
  "SEMA4"
)

target_ec_groups <- c(
  "Adipose Tissue Arterial EC",
  "Adipose Tissue Capillary EC",
  "Adipose Tissue Venous EC"
)

plot_cellchat_pathway_strength <- rankNet(
  cellchat_momac2_ec_merged,
  mode = "comparison",
  comparison = c(1, 2),
  sources.use = "Adipose Tissue Mo-Mac 2",
  targets.use = target_ec_groups,
  signaling = pathways_momac2_to_ec,
  do.stat = TRUE,
  color.use = c("black", "slateblue", "grey")
) +
  aes(x = contribution.scaled, y = name) +
  xlab("Signaling strength") +
  ylab(NULL)

plot_cellchat_pathway_strength

# extract plotting data for source data export and to check pathway p value
df_cellchat_pathway_strength <- plot_cellchat_pathway_strength$data