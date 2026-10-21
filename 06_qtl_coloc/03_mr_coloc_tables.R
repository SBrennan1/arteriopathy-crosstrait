# Joins the MR results to their colocalization results, and compares pQTL and eQTL MR
# directions for gene-outcome pairs significant in both.

source("R/setup.R")
out <- res_dir("qtl_coloc")
fdr <- cfg$qtl_coloc$fdr

pqtl_coloc <- read_tsv(file.path(out, "pqtl_coloc.tsv"), show_col_types = FALSE)
eqtl_coloc <- read_tsv(file.path(out, "eqtl_coloc.tsv"), show_col_types = FALSE)

pqtl <- read_tsv(file.path(cfg$results_dir, "mr/mr_pqtl.tsv"), show_col_types = FALSE) %>%
  left_join(pqtl_coloc %>%
              select(outcome, platform, exposure = protein, chr, start, end, nsnps_coloc = nsnps,
                     abf_pph3, abf_pph4, susie_pph3, susie_pph4, lead_exp, lead_out, coloc_snp,
                     coloc_method = method),
            by = c("outcome", "platform", "exposure")) %>%
  mutate(coloc_evidence = qtl_coloc_evidence(abf_pph4)) %>%
  arrange(outcome, qval, pval)

eqtl <- read_tsv(file.path(cfg$results_dir, "mr/mr_eqtl.tsv"), show_col_types = FALSE) %>%
  left_join(eqtl_coloc %>% select(outcome, tissue, gene_id, nsnps_coloc = nsnps, PP.H3, PP.H4),
            by = c("outcome", "tissue", "gene_id")) %>%
  mutate(coloc_evidence = qtl_coloc_evidence(PP.H4)) %>%
  arrange(outcome, qval, pval)

write_tsv(pqtl, file.path(out, "pqtl_mr_coloc.tsv"))
write_tsv(eqtl, file.path(out, "eqtl_mr_coloc.tsv"))

sig_p <- filter(pqtl, qval < fdr)
sig_e <- filter(eqtl, qval < fdr)
write_counts("qtl_coloc", c(
  pqtl_mr_fdr = nrow(sig_p),
  pqtl_mr_fdr_proteins = n_distinct(sig_p$exposure),
  pqtl_coloc_tested = sum(sig_p$coloc_evidence != "Not tested"),
  pqtl_coloc_strong = sum(sig_p$coloc_evidence == "Strong"),
  pqtl_coloc_suggestive = sum(sig_p$coloc_evidence == "Suggestive"),
  pqtl_coloc_strong_proteins = n_distinct(sig_p$exposure[sig_p$coloc_evidence == "Strong"]),
  eqtl_mr_fdr = nrow(sig_e),
  eqtl_mr_fdr_genes = n_distinct(sig_e$exposure),
  eqtl_coloc_tested = sum(sig_e$coloc_evidence != "Not tested"),
  eqtl_coloc_strong = sum(sig_e$coloc_evidence == "Strong"),
  eqtl_coloc_suggestive = sum(sig_e$coloc_evidence == "Suggestive"),
  eqtl_coloc_strong_genes = n_distinct(sig_e$exposure[sig_e$coloc_evidence == "Strong"])
))

# pQTL vs eQTL direction ----
cross <- sig_p %>%
  select(platform, gene = exposure, outcome, b_pqtl = b, pqtl_coloc = coloc_evidence) %>%
  inner_join(sig_e %>% select(tissue, gene = exposure, outcome, b_eqtl = b, eqtl_coloc = coloc_evidence),
             by = c("gene", "outcome"), relationship = "many-to-many") %>%
  mutate(arterial = tissue %in% cfg$mr$eqtl$arterial_tissues,
         concordant = sign(b_pqtl) == sign(b_eqtl)) %>%
  arrange(gene, outcome, platform, tissue)
write_tsv(cross, file.path(out, "pqtl_eqtl_concordance.tsv"))

cross_sum <- cross %>%
  summarise(comparisons = n(), concordant = sum(concordant),
            pairs = n_distinct(paste(gene, outcome)), .by = c(platform, arterial)) %>%
  bind_rows(cross %>% summarise(comparisons = n(), concordant = sum(concordant),
                                pairs = n_distinct(paste(gene, outcome)), .by = platform)) %>%
  mutate(pct = round(100 * concordant / comparisons, 1))
print(cross_sum)
print(cross %>% filter(!concordant) %>% distinct(gene, outcome, platform) %>% arrange(gene), n = Inf)

write_counts("pqtl_eqtl_concordance", c(
  comparisons = nrow(cross),
  concordant = sum(cross$concordant),
  gene_outcome_pairs = n_distinct(paste(cross$gene, cross$outcome)),
  set_names(cross_sum$comparisons[is.na(cross_sum$arterial)], paste0(cross_sum$platform[is.na(cross_sum$arterial)], "_comparisons")),
  set_names(cross_sum$concordant[is.na(cross_sum$arterial)], paste0(cross_sum$platform[is.na(cross_sum$arterial)], "_concordant"))
))
