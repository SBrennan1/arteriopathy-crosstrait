# Figure 5: (a) expression of target genes across arterial cell types (HCA arterial
# atlas); (b) cell-type specificity relative to locus-matched genes; (c) the same for
# IA target genes in the brain aneurysm atlas.

source("R/setup.R")
library(ggplot2)
library(patchwork)
library(forcats)

p <- cfg$vascular
bd <- p$datasets$brain_aneurysm
q_cut <- 0.05
loci <- read_result("vascular/locus_groups.tsv")
dot_ct <- read_result("vascular/hca_dotplot_celltype.tsv")
enrich <- read_result("vascular/hca_celltype_enrichment.tsv")
brain_spec <- read_result("vascular/brain_ia_specificity.tsv")
brain_donors <- read_result("vascular/brain_donors.tsv")
table1 <- if (file.exists(file.path(cfg$results_dir, "convergence/table1.tsv"))) {
  read_result("convergence/table1.tsv")$Gene
} else character()

theme_set(theme_bw(base_size = 7) +
            theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
                  strip.text = element_text(face = "bold", size = 6.5),
                  legend.key.size = unit(3, "mm"), legend.title = element_text(size = 6.5),
                  plot.tag = element_text(face = "bold", size = 9),
                  plot.title = element_text(size = 7, face = "bold")))
col_up <- "#B2182B"
col_down <- "#2166AC"

hca_lab <- c(`endothelial cell of artery` = "Endothelial cell", `blood vessel smooth muscle cell` = "Smooth muscle cell",
             `microcirculation associated smooth muscle cell` = "Microvascular SMC",
             fibroblast = "Fibroblast", macrophage = "Macrophage", lymphocyte = "Lymphocyte")
brain_lab <- unlist(bd$cell_type_labels)
brain_keep <- setdiff(names(brain_lab), unlist(bd$figure_exclude))
stars <- function(q) case_when(q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ "")

dot_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   axis.text.y = element_text(face = "italic", size = 6),
                   strip.text.y.left = element_text(angle = 0, hjust = 1),
                   strip.placement = "outside", panel.spacing.y = unit(0.8, "mm"),
                   panel.grid.major = element_line(colour = "grey95"))
dot_scales <- list(
  scale_size_area(max_size = 3, labels = scales::percent, name = "Cells expressing"),
  scale_colour_gradient2(low = col_down, mid = "grey92", high = col_up, name = "Scaled expression")
)

# Genes grouped by the cell type of highest scaled expression, ordered within the group
group_genes <- function(d, levels) {
  peak <- d %>%
    slice_max(z, n = 1, with_ties = FALSE, by = gene) %>%
    transmute(gene, group = factor(cell_type, levels), peak_z = z)
  d %>%
    left_join(peak, by = "gene") %>%
    mutate(gene = fct_reorder(gene, peak_z))
}

# Enrichment of target genes over the other genes at the same loci: fold (observed /
# mean of the comparison draws), drawn from 1, with the 95% range of the comparison
# draws as a grey band and BH q-values as stars
enrich_plot <- function(d, y) {
  d %>%
    mutate(fold = observed / null_mean, lo = null_lo / null_mean, hi = null_hi / null_mean,
           call = case_when(q_greater < q_cut ~ "up", q_less < q_cut ~ "down", TRUE ~ "ns"),
           star = stars(pmin(q_greater, q_less))) %>%
    ggplot(aes(y = {{ y }})) +
    geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey85", linewidth = 2.5) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.3) +
    geom_segment(aes(x = 1, xend = fold, yend = {{ y }}, colour = call), linewidth = 0.6) +
    geom_point(aes(x = fold, colour = call), size = 1.8) +
    geom_text(aes(x = Inf, label = star), hjust = 1.15, vjust = 0.75, size = 2.6) +    # stars in a column at the right
    scale_colour_manual(values = c(up = col_up, down = col_down, ns = "grey35"), guide = "none") +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = c(0.05, 0.02))) +
    expand_limits(x = 1 + 1.4 * (max(c(d$observed / d$null_mean, d$null_hi / d$null_mean), na.rm = TRUE) - 1)) +
    labs(x = "Fold enrichment", y = NULL) +
    theme(panel.grid.major.y = element_blank())
}

# a ----
show_genes <- loci %>% filter(source == "ASSET") %>% pull(target_gene) %>% union(table1)
dot <- dot_ct %>%
  filter(gene %in% show_genes, cell_type %in% names(hca_lab)) %>%
  filter(max(pct) >= 0.02, .by = gene) %>%                    # drop near-silent genes
  mutate(z = as.numeric(scale(mean_log)), .by = gene) %>%
  mutate(cell_type = factor(hca_lab[cell_type], levels = unname(hca_lab))) %>%
  group_genes(unname(hca_lab))
pa <- ggplot(dot, aes(cell_type, gene)) +
  geom_point(aes(size = pct, colour = z)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  dot_scales +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL, title = "Arterial atlas") +
  dot_theme +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))

# b ----
pb <- enrich %>%
  filter(cell_type %in% names(hca_lab)) %>%
  mutate(cell_type = factor(hca_lab[cell_type], levels = rev(unname(hca_lab)))) %>%
  enrich_plot(cell_type) +
  labs(title = "Arterial atlas")

# c ----
n_don <- count(brain_donors, region)
pc <- brain_spec %>%
  filter(cell_type %in% brain_keep) %>%
  mutate(condition = str_to_sentence(condition),
         cell_type = factor(brain_lab[cell_type], levels = rev(unname(brain_lab[brain_keep])))) %>%
  enrich_plot(cell_type) +
  facet_wrap(~ condition, nrow = 1) +
  labs(title = "Brain aneurysm atlas") +
  theme(panel.spacing.x = unit(3, "mm"))

n_genes <- n_distinct(dot$gene)
height <- max(130, 45 + n_genes * 3.1)
# The right-hand panels are wrapped so that their axes are not aligned with panel a's
# rotated labels; panel tags are set on each plot for the same reason
right <- wrap_elements(full = (pb + labs(tag = "b")) / (pc + labs(tag = "c")) + plot_layout(heights = c(1, 2.4)))
fig <- ((pa + labs(tag = "a")) | right) + plot_layout(widths = c(1, 1.25))
save_fig(fig, "Figure5.pdf", width = 180, height = height, units = "mm")
