# Shared colours, trait orders and plotting helpers

# Anatomical order (head and neck, intracranial, aorta, coronary, renal/other)
trait_order_anatomical <- c("CeAD", "IA", "AA", "SCAD", "FMD")

# Order and colours used in the MR figures
trait_order_mr <- c("IA", "AA", "FMD", "SCAD", "CeAD")
trait_colors <- c(IA = "#4E79A7", AA = "#F28E2B", FMD = "#59A14F", SCAD = "#E15759", CeAD = "#B07AA1")

status_colors <- c(Known = "#0072B5", Novel = "#BC3C29")
rg_colors <- c(low = "#2166AC", mid = "white", high = "#B2182B")

read_result <- function(path, ...) {
  read_tsv(file.path(cfg$results_dir, path), show_col_types = FALSE, ...)
}

save_fig <- function(plot, file, width, height, ...) {
  ggplot2::ggsave(file.path(res_dir("figures"), file), plot, width = width, height = height, ...)
}

# MR Manhattan panels (Figure 4a/b) ----

mr_block_width <- 1
mr_block_gap <- 0.3

# Best result per gene and outcome, colocalization flag and x positions
mr_manhattan_data <- function(mr, coloc_pairs, fdr = cfg$mr$fdr) {
  mr %>%
    filter(!is.na(pval), pval > 0, trait %in% trait_order_mr) %>%
    group_by(gene, trait) %>%
    slice_min(pval, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(distinct(coloc_pairs, gene, trait) %>% mutate(coloc = TRUE), by = c("gene", "trait")) %>%
    mutate(coloc = coalesce(coloc, FALSE),
           sig = qval < fdr,
           neglog10p = -log10(pval),
           trait = factor(trait, levels = trait_order_mr)) %>%
    arrange(trait, gene) %>%
    group_by(trait) %>%
    mutate(x_pos = (as.numeric(trait) - 1) * (mr_block_width + mr_block_gap) +
             (row_number() - 0.5) / n() * mr_block_width) %>%
    ungroup()
}

mr_manhattan_inputs <- function(strong = cfg$qtl_coloc$strong) {
  pqtl <- read_result("mr/mr_pqtl.tsv") %>%
    transmute(gene = exposure, trait = outcome, platform, beta = b, se, pval, qval)
  pqtl_coloc <- read_result("qtl_coloc/pqtl_coloc.tsv") %>%
    filter(abf_pph4 >= strong | susie_pph4 >= strong) %>%
    transmute(gene = protein, trait = outcome)

  eqtl <- read_result("mr/mr_eqtl.tsv") %>%
    filter(tissue %in% cfg$mr$eqtl$arterial_tissues) %>%
    transmute(gene = exposure, trait = outcome, tissue, beta = b, se, pval, qval)
  eqtl_coloc <- read_result("qtl_coloc/eqtl_coloc.tsv") %>%
    filter(PP.H4 >= strong) %>%
    transmute(gene, trait = outcome)

  list(pqtl = mr_manhattan_data(pqtl, pqtl_coloc), eqtl = mr_manhattan_data(eqtl, eqtl_coloc))
}

# pQTL panel: every colocalising hit plus the top two others per trait.
# eQTL panel: the top five colocalising hits per trait.
mr_labels_pqtl <- function(d, n_extra = 2) {
  bind_rows(
    filter(d, coloc, sig),
    d %>% filter(sig, !coloc) %>% group_by(trait) %>%
      slice_max(neglog10p, n = n_extra, with_ties = FALSE) %>% ungroup()
  ) %>%
    distinct(gene, trait, .keep_all = TRUE) %>%
    filter(!str_detect(gene, "^ENSG"))
}

mr_labels_eqtl <- function(d, n_top = 5) {
  d %>%
    filter(coloc, sig) %>%
    group_by(trait) %>%
    slice_max(neglog10p, n = n_top, with_ties = FALSE) %>%
    ungroup() %>%
    filter(!str_detect(gene, "^ENSG"))
}

# The dashed line is drawn at the largest P value still significant at the FDR threshold
plot_mr_manhattan <- function(data, labels, title) {
  boundary <- data %>% filter(sig) %>% slice_max(pval, n = 1, with_ties = FALSE)
  fdr_line <- if (nrow(boundary) > 0) -log10(boundary$pval) else 3
  max_y <- max(data$neglog10p, na.rm = TRUE) * 1.15

  data <- mutate(data, neglog10p_plot = pmin(neglog10p, max_y * 0.98))
  labels <- mutate(labels, neglog10p_plot = pmin(neglog10p, max_y * 0.98))
  midpoints <- data %>% group_by(trait) %>% summarise(mid = mean(x_pos), .groups = "drop")
  sep_x <- seq_len(length(trait_order_mr) - 1) * (mr_block_width + mr_block_gap) - mr_block_gap / 2

  ggplot(data, aes(x = x_pos, y = neglog10p_plot)) +
    geom_vline(xintercept = sep_x, color = "grey88", linewidth = 0.3) +
    geom_point(data = filter(data, !sig), aes(color = trait),
               size = 1.0, alpha = 0.25, shape = 16) +
    geom_point(data = filter(data, sig, !coloc), aes(color = trait),
               size = 2.3, alpha = 0.75, shape = 21, fill = "white", stroke = 0.7) +
    geom_point(data = filter(data, sig, coloc), aes(fill = trait),
               size = 3.0, shape = 21, color = "grey20", stroke = 0.5) +
    geom_hline(yintercept = fdr_line, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = labels, aes(x = x_pos, y = neglog10p_plot, label = gene),
      size = 3.6, fontface = "italic", max.overlaps = 40,
      segment.size = 0.3, segment.color = "grey50", min.segment.length = 0.1,
      box.padding = 0.5, point.padding = 0.3, force = 3, force_pull = 0.3, seed = 42
    ) +
    scale_x_continuous(breaks = midpoints$mid, labels = midpoints$trait,
                       limits = c(-mr_block_gap / 2, max(data$x_pos) + mr_block_gap / 2),
                       expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)), limits = c(0, max_y)) +
    scale_color_manual(values = trait_colors, guide = "none") +
    scale_fill_manual(values = trait_colors, guide = "none") +
    labs(x = NULL, y = expression(-log[10](italic(P))), title = title) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line.y = element_line(colour = "black", linewidth = 0.4),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.4),
      axis.text.x = element_text(size = 12, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 11, colour = "black"),
      axis.title.y = element_text(size = 12, colour = "black"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
}
