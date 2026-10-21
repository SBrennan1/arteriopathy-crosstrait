# Figure S3: (a) genetic correlation |Z| for arteriopathy pairs and matched controls;
# (b) strong colocalization at the ASSET loci for arteriopathy pairs, matched controls
# and hair color.

source("R/setup.R")
library(ggplot2)
library(patchwork)

rg_ctrl <- read_result("negative_controls/rg_controls.tsv")
rg_pairs <- read_result("negative_controls/rg_vs_controls.tsv")
coloc_draw <- read_result("negative_controls/coloc_null_by_draw.tsv")
coloc_ref <- read_result("negative_controls/coloc_reference.tsv")
coloc_obs <- read_result("coloc/cross_trait_coloc.tsv")
draws <- fread(cfg$negative_controls$draws)

lab <- function(x) str_replace_all(x, c("^AA$" = "AAD"))
pair_lab <- function(a, b) paste(lab(a), lab(b), sep = "-")
cols <- c(Arteriopathy = "#BC3C29", `Matched controls` = "grey60", `Hair color` = "#2E6FB7")

theme_set(theme_bw(base_size = 7) +
            theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                  plot.tag = element_text(face = "bold", size = 9),
                  legend.key.size = unit(3, "mm"), legend.title = element_blank(),
                  legend.position = "bottom"))

# a ----
matched <- distinct(draws[ldsc_subset == TRUE, .(matched_to = trait, id)])
rg_null <- rg_pairs %>%
  select(trait1, trait2) %>%
  mutate(null = map2(trait1, trait2, function(a, b) {
    bind_rows(rg_ctrl %>% filter(trait == a, id %in% matched[matched_to == b]$id),
              rg_ctrl %>% filter(trait == b, id %in% matched[matched_to == a]$id))
  })) %>%
  unnest(null) %>%
  mutate(pair = pair_lab(trait1, trait2), abs_z = abs(rg / se), group = "Matched controls")
obs <- rg_pairs %>% mutate(pair = pair_lab(trait1, trait2), abs_z = abs(rg / se), group = "Arteriopathy")
pair_order <- obs %>% arrange(desc(abs_z)) %>% pull(pair)

pa <- ggplot(mapping = aes(factor(pair, levels = pair_order), abs_z)) +
  geom_hline(yintercept = qnorm(0.975), linetype = 2, colour = "grey50", linewidth = 0.3) +
  geom_boxplot(data = rg_null, aes(colour = group), outlier.size = 0.3, width = 0.6, linewidth = 0.3) +
  geom_point(data = obs, aes(colour = group), size = 2) +
  scale_colour_manual(values = cols, breaks = c("Arteriopathy", "Matched controls"),
                      labels = c("Arteriopathy pair", "Matched controls")) +
  labs(x = NULL, y = "Genetic correlation |Z|") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# b ----
rates <- tibble(
  group = factor(names(cols), levels = names(cols)),
  tests = c(sum(!is.na(coloc_obs$PP.H4.abf)), sum(coloc_draw$pairs_tested), sum(!is.na(coloc_ref$PP.H4.abf))),
  strong = c(sum(coloc_obs$PP.H4.abf >= 0.8, na.rm = TRUE), sum(coloc_draw$strong_pairs),
             sum(coloc_ref$PP.H4.abf >= 0.8, na.rm = TRUE))
) %>%
  mutate(pct = 100 * strong / tests,
         label = paste0(prettyNum(strong, big.mark = ","), "/", prettyNum(tests, big.mark = ",")))

pb <- ggplot(rates, aes(group, pct, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label), vjust = -0.4, size = 2.3) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Tests with strong colocalization (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p <- (pa + labs(tag = "a") | pb + labs(tag = "b")) + plot_layout(widths = c(3, 1.3))
save_fig(p, "FigureS3.pdf", width = 180, height = 85, units = "mm")
