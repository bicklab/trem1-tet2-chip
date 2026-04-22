library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(purrr)

# =========================
# Figure: TREM1+ fraction by cell type and CHIP status
# =========================

# ---- user-defined settings ----
target_gene <- "TREM1"

celltype_order <- c(
  "Adipose Tissue cDC2B",
  "Adipose Tissue IM",
  "Adipose Tissue LAM",
  "Adipose Tissue Mo-Mac 2"
)

chip_order <- c("Control", "TET2")

chip_colors <- c(
  "Control" = "#6baed6",
  "TET2"    = "#e34a33"
)

# ---- subset Seurat object to relevant cell types ----
seu_trem1_fraction <- subset(
  obj,
  subset = Cell_Type %in% celltype_order
)

# ---- extract expression and metadata ----
trem1_expr_vec <- GetAssayData(
  seu_trem1_fraction,
  assay = "RNA",
  slot = "data"
)[target_gene, ]

df_trem1_cell_level <- seu_trem1_fraction@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, Cell_Type, CHIP) %>%
  mutate(
    trem1_positive = trem1_expr_vec > 0,
    Cell_Type = factor(Cell_Type, levels = celltype_order),
    CHIP = factor(CHIP, levels = chip_order)
  )

# ---- summarise counts and fractions ----
df_trem1_positive_fraction <- df_trem1_cell_level %>%
  group_by(Cell_Type, CHIP, trem1_positive) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(
    Cell_Type,
    CHIP,
    trem1_positive,
    fill = list(n = 0)
  ) %>%
  group_by(Cell_Type, CHIP) %>%
  mutate(
    total_cells = sum(n),
    pct = 100 * n / total_cells
  ) %>%
  ungroup() %>%
  filter(trem1_positive) %>%
  mutate(
    n_label = paste0("n=", n),
    pct_label = sprintf("%.1f%%", pct)
  )

# ---- statistical testing: Control vs TET2 within each cell type ----
# 2x2 Fisher's exact test on TREM1+ vs TREM1- by CHIP
df_trem1_positive_stats <- df_trem1_cell_level %>%
  group_by(Cell_Type) %>%
  group_modify(~{
    tab <- table(.x$CHIP, .x$trem1_positive)
    
    # force both dimensions to exist
    all_rows <- chip_order
    all_cols <- c(FALSE, TRUE)
    
    mat <- matrix(
      0,
      nrow = length(all_rows),
      ncol = length(all_cols),
      dimnames = list(all_rows, all_cols)
    )
    
    existing_rows <- intersect(rownames(tab), all_rows)
    existing_cols <- intersect(colnames(tab), c("FALSE", "TRUE"))
    mat[existing_rows, existing_cols] <- tab[existing_rows, existing_cols]
    
    fisher_res <- fisher.test(mat)
    
    tibble(
      p = fisher_res$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    stars = case_when(
      p_adj <= 0.0001 ~ "****",
      p_adj <= 0.001  ~ "***",
      p_adj <= 0.01   ~ "**",
      p_adj <= 0.05   ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  left_join(
    df_trem1_positive_fraction %>%
      group_by(Cell_Type) %>%
      summarise(y_max = max(pct), .groups = "drop"),
    by = "Cell_Type"
  ) %>%
  mutate(
    xmin = 1,
    xmax = 2,
    y.position = y_max + 8
  )

# ---- make plotting-friendly facet-specific x positions ----
df_trem1_positive_fraction <- df_trem1_positive_fraction %>%
  mutate(
    CHIP = factor(CHIP, levels = chip_order)
  )

# ---- plot ----
plot_trem1_positive_fraction <- ggplot(
  df_trem1_positive_fraction,
  aes(x = CHIP, y = pct, fill = CHIP)
) +
  geom_col(
    color = "black",
    width = 0.7
  ) +
  geom_text(
    aes(label = n_label),
    vjust = 1.3,
    color = "white",
    size = 4
  ) +
  geom_text(
    aes(label = pct_label),
    vjust = -0.5,
    size = 4
  ) +
  facet_wrap(~ Cell_Type, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = chip_colors) +
  scale_y_continuous(
    name = "TREM1+ cells (% of total)",
    expand = expansion(mult = c(0, 0.18))
  ) +
  xlab(NULL) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    strip.placement = "outside",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggpubr::stat_pvalue_manual(
    df_trem1_positive_stats,
    label = "stars",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.01,
    size = 5
  )

plot_trem1_positive_fraction
