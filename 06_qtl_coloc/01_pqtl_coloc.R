# Colocalization of cis-pQTLs with the outcome GWAS for each FDR-significant pQTL MR
# association (gene +/- 500 kb): coloc.abf, with coloc.susie where 1000G EUR LD allows.

source("R/setup.R")
p <- cfg$qtl_coloc
out <- res_dir("qtl_coloc")

jobs <- read_tsv(file.path(cfg$results_dir, "mr/mr_pqtl.tsv"), show_col_types = FALSE) %>%
  filter(qval < p$fdr) %>%
  distinct(outcome, platform, protein = exposure)
jobs <- jobs %>% left_join(protein_regions(unique(jobs$protein)), by = "protein")

# Outcome GWAS restricted to the tested windows
gwas_cols <- c("rsid", "chr", "pos", "effect_allele", "other_allele", "eaf", "beta", "se", "pval")
gwas <- map(set_names(unique(jobs$outcome)), function(t) {
  win <- jobs %>% filter(outcome == t, !is.na(chr)) %>% distinct(chr, start, end)
  read_univar(t, gwas_cols) %>%
    rows_in_windows(win) %>%
    drop_na(rsid, beta, se, effect_allele, other_allele) %>%
    mutate(pval = pmax(pval, 1e-300))
})

read_exposure <- function(platform, protein, region) {
  switch(platform,
         Olink = read_olink_region(protein, region),
         deCODE = read_decode_region(protein, region),
         stop("unknown platform ", platform))
}

# Run ----
res <- pmap_dfr(jobs, function(outcome, platform, protein, chr, start, end) {
  r <- tryCatch({
    if (is.na(chr)) stop("no gene coordinates")
    region <- list(chr = chr, start = start, end = end)
    outcome_region <- gwas[[outcome]] %>%
      filter(chr == region$chr, pos >= region$start, pos <= region$end) %>%
      distinct(rsid, .keep_all = TRUE)
    coloc_pqtl(read_exposure(platform, protein, region), outcome_region, region,
               n_exp = p[[tolower(platform)]]$sample_size,
               ncase = cfg$traits[[outcome]]$ncase, ncontrol = cfg$traits[[outcome]]$ncontrol) %>%
      mutate(reason = NA_character_)
  }, error = function(e) tibble(reason = conditionMessage(e)))
  bind_cols(tibble(outcome, platform, protein), r)
})

write_tsv(res %>% filter(!is.na(reason)) %>% select(outcome, platform, protein, reason),
          file.path(out, "pqtl_coloc_not_run.tsv"))

write_tsv(res %>%
            filter(is.na(reason)) %>%
            select(outcome, platform, protein, chr, start, end, nsnps, susie_pph3, susie_pph4,
                   abf_pph3, abf_pph4, lead_exp, lead_out, coloc_snp, method) %>%
            arrange(outcome, platform, protein),
          file.path(out, "pqtl_coloc.tsv"))
