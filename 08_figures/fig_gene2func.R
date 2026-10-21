# Figure S4: FUMA GENE2FUNC tissue expression (GTEx v8) of the ASSET genes and
# enrichment among tissue DEG sets.

source("R/setup.R")
library(ggplot2)
library(cowplot)
library(grid)
library(ComplexHeatmap)

g2f <- cfg$enrichment$gene2func_dir

# Expression heatmap ----
exp_raw <- read.delim(file.path(g2f, "gtex_v8_ts_avg_normTPM_exp.txt"), check.names = FALSE)
exp_mat <- as.matrix(exp_raw[, setdiff(names(exp_raw), c("ensg", "symbol"))])
rownames(exp_mat) <- make.unique(exp_raw$symbol)

# z-score per gene across tissues, capped at +/- 3
exp_z <- t(scale(t(exp_mat)))
colnames(exp_z) <- gsub("_", " ", gsub("^.*ts_general_", "", colnames(exp_z)))
exp_z[exp_z > 3] <- 3
exp_z[exp_z < -3] <- -3

ht <- Heatmap(
  exp_z,
  name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), c("#4575B4", "#FFFFBF", "#D73027")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_dend_width = unit(15, "mm"),
  column_dend_height = unit(15, "mm"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(25, "mm"),
    legend_width = unit(4, "mm")
  ),
  border = TRUE,
  width = unit(160, "mm"),
  height = unit(240, "mm")
)

ht_grob <- grid.grabExpr(
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       padding = unit(c(5, 5, 5, 5), "mm")),
  width = 8, height = 12
)

# DEG enrichment ----
# Bar length: -log10 P; colour: adjusted P
deg <- read.delim(file.path(g2f, "gtex_v8_ts_DEG.txt"), check.names = FALSE) %>%
  mutate(
    tissue = gsub("_", " ", GeneSet),
    direction = case_when(
      grepl("up", Category, ignore.case = TRUE) ~ "Up-regulated DEG",
      grepl("down", Category, ignore.case = TRUE) ~ "Down-regulated DEG",
      TRUE ~ "Other"
    ),
    neglog10p = -log10(p),
    bar_fill = if_else(adjP < 0.05, "#D73027", "#7BAFD4")
  ) %>%
  filter(direction != "Other")

tissue_order <- deg %>%
  group_by(tissue) %>%
  summarise(max_enrich = max(neglog10p), .groups = "drop") %>%
  arrange(max_enrich) %>%
  pull(tissue)
deg <- mutate(deg, tissue = factor(tissue, levels = tissue_order))

p_deg <- ggplot(deg, aes(x = neglog10p, y = tissue)) +
  geom_col(aes(fill = bar_fill), width = 0.7) +
  scale_fill_identity() +
  facet_wrap(~ direction, ncol = 1, scales = "free_x") +
  labs(x = expression(-log[10](italic(P))), y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 6, colour = "grey30"),
    strip.text = element_text(size = 7, face = "bold", hjust = 0),
    strip.background = element_blank(),
    panel.spacing = unit(8, "mm"),
    plot.margin = margin(35, 10, 5, 5)
  )

fig <- plot_grid(ht_grob, p_deg, labels = c("A", "B"), label_size = 12, label_fontface = "bold",
                 ncol = 2, rel_widths = c(1.6, 1), align = "none")
save_fig(fig, "FigureS4.pdf", width = 14, height = 13)
