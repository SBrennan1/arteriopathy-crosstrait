# Figure S2: regional plots (+/- 500 kb) of the novel MTAG loci. Requires LDLINK_TOKEN.

source("R/setup.R")
library(ggplot2)
library(cowplot)
library(locuszoomr)
select <- dplyr::select   # masked by AnnotationDbi via locuszoomr

token <- Sys.getenv("LDLINK_TOKEN")
if (!nzchar(token)) stop("Set LDLINK_TOKEN to an LDlink API token")

novel <- read_result("mtag/mtag_leads.tsv") %>%
  filter(novel) %>%
  left_join(read_result("v2g/mtag_v2g.tsv") %>% select(rsid, gene), by = "rsid") %>%
  arrange(target_trait, p_mtag) %>%
  mutate(panel_label = paste0(gene, " (", target_trait, ")"))

# MTAG summary statistics with the column names locuszoomr expects
mtag <- map(set_names(unique(novel$target_trait)), function(t) {
  vroom(file.path(cfg$results_dir, "mtag/fuma_input", paste0(t, "_mtag.txt.gz")),
        col_select = c(rsid, chr, pos, pval), show_col_types = FALSE) %>%
    transmute(rsid, chrom = as.integer(chr), pos = as.integer(pos), p = as.numeric(pval),
              logp = -log10(p)) %>%
    as.data.frame()
})

ens_db <- ensembldb::addFilter(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                               AnnotationFilter::TxBiotypeFilter("protein_coding"))

plots <- map(seq_len(nrow(novel)), function(i) {
  l <- novel[i, ]
  loc <- locus(data = mtag[[l$target_trait]], index_snp = l$rsid, flank = 500e3, ens_db = ens_db)
  loc <- link_LD(loc, token = token, pop = "EUR")
  locus_ggplot(loc, labels = "index", pcutoff = 5e-8, legend_pos = "topright", maxrows = 3) +
    ggtitle(l$panel_label) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0),
          text = element_text(size = 8),
          plot.margin = margin(t = 14, r = 8, b = 6, l = 14, unit = "pt"))
})

n_col <- min(3, length(plots))
n_row <- ceiling(length(plots) / 3)
fig <- plot_grid(plotlist = plots, ncol = n_col, nrow = n_row,
                 labels = LETTERS[seq_along(plots)], label_size = 12, label_fontface = "bold",
                 label_x = 0, label_y = 1, hjust = -0.2, vjust = 1.2, align = "hv", axis = "tblr")

save_fig(fig, "FigureS2.pdf", width = 100 * n_col, height = 110 * n_row, units = "mm", dpi = 600)
