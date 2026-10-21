# Figure 4: proteome-wide (a) and arterial transcriptome-wide (b) MR. Filled points
# colocalize (PP.H4 >= 0.8).

source("R/setup.R")
library(ggplot2)
library(ggrepel)
library(cowplot)

mh <- mr_manhattan_inputs()

# Manhattan panels ----
panel_a <- plot_mr_manhattan(mh$pqtl, mr_labels_pqtl(mh$pqtl), "a")
panel_b <- plot_mr_manhattan(mh$eqtl, mr_labels_eqtl(mh$eqtl), "b")

fig <- plot_grid(panel_a, panel_b, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1))
save_fig(fig, "Figure4.pdf", width = 7.2, height = 8.5, device = cairo_pdf)
