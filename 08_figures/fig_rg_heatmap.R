# Figure 2: (a) pairwise genetic correlation (LDSC; * BH FDR < 0.05) and (b) pairwise
# colocalization at the ASSET loci (strong / tested).

source("R/setup.R")
library(ggplot2)
library(patchwork)

lab <- function(x) recode(x, AA = "AAD")
lvls <- lab(trait_order_anatomical)

# Lower triangle in anatomical order
lower_tri <- function(d) {
  bind_rows(d, rename(d, trait1 = trait2, trait2 = trait1)) %>%
    mutate(trait1 = factor(lab(trait1), levels = lvls), trait2 = factor(lab(trait2), levels = lvls)) %>%
    filter(as.numeric(trait1) > as.numeric(trait2))
}

heat_theme <- theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 8, colour = "black"), panel.grid = element_blank(),
        legend.key.height = unit(4, "mm"), legend.key.width = unit(2.5, "mm"),
        legend.title = element_text(size = 7.5), legend.text = element_text(size = 7),
        plot.tag = element_text(face = "bold", size = 10))

# a ----
rg <- read_result("ldsc/rg_pairwise.tsv") %>%
  transmute(trait1 = std_trait(p1), trait2 = std_trait(p2), rg, se, p) %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH"))
rg_df <- lower_tri(rg) %>%
  mutate(label = paste0(sprintf("%.2f", rg), if_else(coalesce(p_adj_BH < 0.05, FALSE), "*", "")))

pa <- ggplot(rg_df, aes(trait2, trait1, fill = rg)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.8) +
  scale_fill_gradient2(name = expression(r[g]), limits = c(-1, 1), midpoint = 0,
                       low = rg_colors[["low"]], mid = rg_colors[["mid"]], high = rg_colors[["high"]]) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  heat_theme

# b ----
coloc_df <- read_result("coloc/cross_trait_coloc.tsv") %>%
  filter(!is.na(PP.H4.abf)) %>%
  mutate(trait1 = std_trait(trait1), trait2 = std_trait(trait2),
         a = pmin(trait1, trait2), b = pmax(trait1, trait2)) %>%
  summarise(tested = n(), strong = sum(PP.H4.abf >= 0.8), .by = c(a, b)) %>%
  rename(trait1 = a, trait2 = b)
all_pairs <- as_tibble(t(combn(sort(trait_order_anatomical), 2)), .name_repair = ~ c("trait1", "trait2"))
coloc_df <- all_pairs %>%
  left_join(coloc_df, by = c("trait1", "trait2")) %>%
  mutate(tested = coalesce(tested, 0L), strong = coalesce(strong, 0L),
         share = if_else(tested > 0, strong / tested, NA_real_)) %>%
  lower_tri() %>%
  mutate(label = paste0(strong, "/", tested))

pb <- ggplot(coloc_df, aes(trait2, trait1, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.8) +
  scale_fill_gradient(name = "Share of loci\ncolocalized", low = "#F4F4F4", high = "#3B6E8F",
                      limits = c(0, 1), labels = scales::percent, na.value = "grey90") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  heat_theme

fig <- (pa + labs(tag = "a")) + (pb + labs(tag = "b")) + plot_layout(nrow = 1)
save_fig(fig, "Figure2.pdf", width = 180, height = 75, units = "mm", bg = "white")

print(coloc_df %>% select(trait1, trait2, strong, tested))
