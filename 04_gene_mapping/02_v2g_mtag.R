# Variant-to-gene mapping of the pairwise MTAG lead SNPs.

source("R/setup.R")
out <- res_dir("v2g")

leads <- read_tsv(file.path(cfg$results_dir, "mtag/mtag_leads_all.tsv"), show_col_types = FALSE) %>%
  group_by(rsid) %>%
  summarise(chr = first(chr), pos = first(pos),
            a1 = first(na.omit(effect_allele_univ)), a2 = first(na.omit(other_allele_univ)),
            traits = paste(unique(target_trait), collapse = ","), .groups = "drop")

res <- v2g_leads(leads, "MTAG")
v2g <- res$v2g %>% arrange(chr, pos)

write_tsv(v2g, file.path(out, "mtag_v2g.tsv"))
write_tsv(res$coloc, file.path(out, "mtag_gtex_coloc.tsv"))

ev <- count(v2g, evidence)
write_counts("v2g_mtag", c(
  lead_snps = nrow(v2g),
  setNames(ev$n, paste0("evidence_", str_replace_all(tolower(ev$evidence), "[^a-z0-9]+", "_"))),
  functional_support = sum(v2g$evidence != "Nearest gene"),
  overridden = sum(v2g$overridden)
))
