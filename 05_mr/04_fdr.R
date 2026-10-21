# Storey q-values within each platform (pQTL) or tissue (eQTL) and outcome
# (Benjamini-Hochberg where pi0 cannot be estimated).

source("R/setup.R")
library(qvalue)
out <- res_dir("mr")
fdr <- cfg$mr$fdr

apply_qvalue <- function(d, groups) {
  d %>%
    group_by(across(all_of(groups))) %>%
    mutate(qval = tryCatch(qvalue(pval)$qvalues,
                           error = function(e) p.adjust(pval, method = "BH"))) %>%
    ungroup()
}

read_mr <- function(file) {
  read_tsv(file.path(out, file), show_col_types = FALSE) %>%
    mutate(outcome = std_trait(outcome)) %>%
    filter(outcome %in% traits)
}

pqtl <- read_mr("mr_pqtl_unadjusted.tsv") %>%
  apply_qvalue(c("platform", "outcome")) %>%
  select(platform, exposure, outcome, nsnp, method, b, se, pval, qval, snps) %>%
  arrange(platform, outcome, pval)

eqtl <- read_mr("mr_eqtl_unadjusted.tsv") %>%
  apply_qvalue(c("tissue", "outcome")) %>%
  select(tissue, exposure, gene_id, outcome, nsnp, method, b, se, pval, qval, snps) %>%
  arrange(tissue, outcome, pval)

write_tsv(pqtl, file.path(out, "mr_pqtl.tsv"))
write_tsv(eqtl, file.path(out, "mr_eqtl.tsv"))

# Counts ----
tally <- function(d, prefix, by) {
  s <- d %>%
    group_by(across(all_of(by))) %>%
    summarise(tested = n(), fdr05 = sum(qval < fdr, na.rm = TRUE), .groups = "drop") %>%
    unite("key", all_of(by)) %>%
    pivot_longer(c(tested, fdr05))
  set_names(s$value, paste(prefix, s$key, s$name, sep = "_"))
}

arterial <- filter(eqtl, tissue %in% cfg$mr$eqtl$arterial_tissues)
sig_proteins <- pqtl %>% filter(qval < fdr) %>% group_by(platform) %>%
  summarise(n = n_distinct(exposure))

write_counts("mr", c(
  tally(pqtl, "pqtl", c("platform", "outcome")),
  set_names(sig_proteins$n, paste0("pqtl_", sig_proteins$platform, "_proteins_fdr05")),
  tally(eqtl, "eqtl", "outcome"),
  tally(arterial, "eqtl_arterial", "outcome"),
  eqtl_arterial_genes_fdr05 = n_distinct(arterial$exposure[which(arterial$qval < fdr)])
))
