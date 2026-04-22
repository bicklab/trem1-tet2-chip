library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tibble)
library(readr)

# ============================================================
# Figure: TREM1 expression in neutrophils by condition
# Statistical test: Wilcoxon rank-sum test
# ============================================================

# ----------------------------
# paths
# ----------------------------
input_rds <- "/home/rstudio/trem1/neutrophil_v2_012026/rds_object/2_merged_neut_annotated_obj.rds"
fig_dir   <- "/home/rstudio/trem1/neutrophil_v2_012026/figures"
data_dir  <- "/home/rstudio/trem1/neutrophil_v2_012026/source_data"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# load object
# ----------------------------
seu_neutrophils <- readRDS(input_rds)

# ----------------------------
# extract plotting data
# ----------------------------
df_trem1 <- FetchData(
  seu_neutrophils,
  vars = c("TREM1", "condition")
) %>%
  as_tibble() %>%
  filter(!is.na(condition)) %>%
  mutate(
    condition = factor(condition, levels = c("Control", "TET2"))
  )

# ----------------------------
# Wilcoxon test
# ----------------------------
wilcox_res <- wilcox.test(TREM1 ~ condition, data = df_trem1)

df_trem1_stats <- tibble(
  gene = "TREM1",
  group1 = "Control",
  group2 = "TET2",
  p_value = wilcox_res$p.value
)

write_tsv(df_trem1, file.path(data_dir, "panelE_trem1_violin_source_data.tsv"))
write_tsv(df_trem1_stats, file.path(data_dir, "panelE_trem1_violin_stats.tsv"))

# ----------------------------
# plot
# ----------------------------
plot_panelE_trem1 <- ggplot(
  df_trem1,
  aes(x = condition, y = TREM1, fill = condition)
) +
  geom_violin(trim = FALSE, color = "black", width = 0.9) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.25) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  scale_fill_manual(values = c("Control" = "#F8766D", "TET2" = "#00BFC4")) +
  labs(
    x = NULL,
    y = "Expression level",
    title = "TREM1"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0)
  )

ggsave(
  filename = file.path(fig_dir, "panelE_trem1_violin.pdf"),
  plot = plot_panelE_trem1,
  width = 4.2,
  height = 4.0
)

ggsave(
  filename = file.path(fig_dir, "panelE_trem1_violin.png"),
  plot = plot_panelE_trem1,
  width = 4.2,
  height = 4.0,
  dpi = 300
)