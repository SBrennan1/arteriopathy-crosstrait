# Figure S1: pairwise MTAG (top) against the univariate GWAS (bottom) for each MTAG
# target trait.

source("R/setup.R")
library(ggplot2)
library(patchwork)
library(ggrepel)

gws <- 5e-8

# GRCh37 chromosome lengths
chr_info <- tibble(
  chr = 1:22,
  chr_len = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
              146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
              102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566)
) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0), chr_center = chr_start + chr_len / 2)
x_max <- max(chr_info$chr_start + chr_info$chr_len)

leads <- read_result("mtag/mtag_leads.tsv") %>%
  left_join(read_result("v2g/mtag_v2g.tsv") %>% select(rsid, gene), by = "rsid") %>%
  select(rsid, chr, pos, target_trait, gene, p_mtag, p_univ, novel_mtag)

# All SNPs with P < p_keep, every thin-th SNP (by position) of the rest
prep_gwas <- function(d, p_keep = 1e-4, thin = 40) {
  d <- d %>%
    filter(!is.na(chr), !is.na(pos), !is.na(pval), chr %in% 1:22) %>%
    mutate(chr = as.integer(chr), neglog10p = -log10(pval)) %>%
    left_join(select(chr_info, chr, chr_start), by = "chr") %>%
    mutate(cumpos = pos + chr_start)
  bind_rows(
    filter(d, pval < p_keep),
    d %>% filter(pval >= p_keep) %>% group_by(chr) %>% arrange(pos) %>%
      slice(seq(1, n(), by = thin)) %>% ungroup()
  )
}

x_scale <- scale_x_continuous(breaks = chr_info$chr_center, labels = chr_info$chr,
                              limits = c(0, x_max), expand = c(0.01, 0))
# Alternating chromosome colours and the label colours share one colour scale
chr_colors <- c(`0` = "#1F77B4", `1` = "#7FAFDB")

miami_plot <- function(trait, raster_dpi = 150) {
  mtag <- vroom(file.path(cfg$results_dir, "mtag/fuma_input", paste0(trait, "_mtag.txt.gz")),
                col_select = c(rsid, chr, pos, pval), show_col_types = FALSE) %>%
    prep_gwas()
  univ <- read_univar(trait, cols = c("rsid", "chr", "pos", "pval")) %>%
    prep_gwas()

  lead <- leads %>%
    filter(target_trait == trait) %>%
    left_join(select(chr_info, chr, chr_start), by = "chr") %>%
    mutate(cumpos = pos + chr_start,
           point_type = if_else(novel_mtag, "Novel", "Known"))
  lead_mtag <- mutate(lead, neglog10p = -log10(p_mtag))
  lead_univ <- mutate(lead, neglog10p = -log10(p_univ))

  # One label per gene, at its most significant lead
  to_label <- lead_mtag %>%
    filter(!grepl("^RP[0-9]|^AC[0-9]|^AP[0-9]|^LINC|^LOC[0-9]|^RPS|^RPL", gene)) %>%
    group_by(gene) %>%
    slice_min(p_mtag, n = 1, with_ties = FALSE) %>%
    ungroup()

  max_y_mtag <- min(max(c(mtag$neglog10p, lead_mtag$neglog10p), na.rm = TRUE) * 1.15, 50)
  max_y_univ <- min(max(c(univ$neglog10p, lead_univ$neglog10p), na.rm = TRUE) * 1.15, 50)

  p_mtag <- ggplot() +
    ggrastr::rasterise(
      geom_point(data = mtag, aes(x = cumpos, y = neglog10p, color = factor(chr %% 2)),
                 size = 0.3, alpha = 0.6),
      dpi = raster_dpi) +
    geom_hline(yintercept = -log10(gws), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_point(data = lead_mtag, aes(x = cumpos, y = neglog10p, fill = point_type),
               color = "black", size = 2.5, shape = 21, stroke = 0.5) +
    geom_text_repel(data = to_label,
                    aes(x = cumpos, y = neglog10p, label = gene, color = point_type),
                    size = 2.5, fontface = "italic", max.overlaps = 20,
                    nudge_y = max_y_mtag * 0.05, segment.size = 0.2, segment.color = "grey50",
                    min.segment.length = 0, box.padding = 0.3, force = 2, show.legend = FALSE) +
    scale_color_manual(values = c(chr_colors, status_colors), guide = "none") +
    scale_fill_manual(values = status_colors, name = NULL) +
    x_scale +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, max_y_mtag)) +
    labs(y = expression(-log[10](P) ~ "MTAG"), x = NULL, title = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 10, 0, 10),
      legend.position = "none"
    )

  p_univ <- ggplot() +
    ggrastr::rasterise(
      geom_point(data = univ, aes(x = cumpos, y = neglog10p, color = factor(chr %% 2)),
                 size = 0.3, alpha = 0.6),
      dpi = raster_dpi) +
    scale_color_manual(values = chr_colors, guide = "none") +
    geom_hline(yintercept = -log10(gws), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_point(data = filter(lead_univ, p_univ < gws), aes(x = cumpos, y = neglog10p),
               fill = status_colors[["Known"]], color = "black", size = 2, shape = 21, stroke = 0.5) +
    x_scale +
    scale_y_reverse(expand = expansion(mult = c(0.05, 0)), limits = c(max_y_univ, 0)) +
    labs(y = expression(-log[10](P) ~ "Univariate"), x = "Chromosome") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(size = 7),
      panel.grid = element_blank(),
      plot.margin = margin(0, 10, 5, 10)
    )

  p_mtag / p_univ +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(title = trait,
                    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)))
}

miami <- list()
for (t in names(cfg$mtag$outputs)) {
  miami[[t]] <- miami_plot(t)
}

panel <- map(miami, wrap_elements)
combined <- (panel$CeAD | panel$SCAD) / (panel$IA | panel$AA) +
  plot_annotation(title = "Pairwise MTAG vs Univariate GWAS",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
save_fig(combined, "FigureS1.pdf", width = 18, height = 10, dpi = 300)
