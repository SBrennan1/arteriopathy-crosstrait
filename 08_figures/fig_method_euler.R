# Figure S5: overlap of loci identified by ASSET, CPASSOC and pairwise MTAG, labelled
# by mapped gene.

source("R/setup.R")
library(eulerr)

loci <- read_result("overlap/locus_by_method.tsv")

v2g <- bind_rows(
  read_result("v2g/asset_v2g.tsv") %>% transmute(rsid = locus_rsid, gene, priority = 1),
  read_result("v2g/cpassoc_v2g.tsv") %>% transmute(rsid, gene, priority = 2),
  read_result("v2g/mtag_v2g.tsv") %>% transmute(rsid, gene, priority = 3)
) %>%
  filter(!is.na(gene))

# Within a method, the locus lead variant's gene is preferred
locus_gene <- loci %>%
  select(locus, lead_rsid, rsid = rsids) %>%
  separate_rows(rsid, sep = ",") %>%
  inner_join(v2g, by = "rsid") %>%
  arrange(locus, priority, rsid != lead_rsid) %>%
  distinct(locus, .keep_all = TRUE) %>%
  select(locus, gene)

loci <- loci %>%
  left_join(locus_gene, by = "locus") %>%
  mutate(gene = coalesce(gene, lead_rsid))

sets <- list(
  ASSET = unique(loci$gene[loci$ASSET]),
  CPASSOC = unique(loci$gene[loci$CPASSOC]),
  MTAG = unique(loci$gene[loci$MTAG])
)
fit <- euler(sets)

pdf(file.path(res_dir("figures"), "FigureS5.pdf"), width = 8, height = 6)
print(plot(fit,
           fills = list(fill = c("#d73027", "#4575b4", "#1a9850"), alpha = 0.5),
           labels = list(fontsize = 14),
           quantities = list(fontsize = 12)))
invisible(dev.off())
