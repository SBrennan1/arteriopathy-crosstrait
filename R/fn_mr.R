# Two-sample MR (TwoSampleMR): instruments matched by rsID, Wald ratio for one SNP and
# fixed-effect IVW for two or more.

read_mr_outcome <- function(trait, rsids = NULL) {
  d <- read_univar(trait, cols = c("rsid", "effect_allele", "other_allele", "eaf", "beta", "se", "pval")) %>%
    mutate(rsid = na_if(trimws(as.character(rsid)), ""))
  if (!is.null(rsids)) d <- filter(d, rsid %in% rsids)
  d
}

# ivs: one row per instrument (exposure, rsid, beta, se, effect_allele, other_allele,
# pval, optionally eaf)
run_mr <- function(ivs, outcome, outcome_name) {
  exp <- ivs %>%
    filter(!is.na(rsid), rsid %in% outcome$rsid) %>%
    distinct(exposure, rsid, .keep_all = TRUE) %>%
    select(exposure, SNP = rsid, beta, se, effect_allele, other_allele, any_of("eaf"), pval)
  if (nrow(exp) == 0) return(NULL)

  exp_dat <- TwoSampleMR::format_data(as.data.frame(exp), type = "exposure",
                                      phenotype_col = "exposure", snp_col = "SNP")
  exp_dat$id.exposure <- exp_dat$exposure

  out <- outcome %>% filter(rsid %in% exp_dat$SNP) %>% mutate(outcome = outcome_name)
  out_dat <- TwoSampleMR::format_data(as.data.frame(out), type = "outcome",
                                      phenotype_col = "outcome", snp_col = "rsid")
  out_dat$id.outcome <- outcome_name

  harmonise <- function(e) TwoSampleMR::harmonise_data(e, out_dat, action = cfg$mr$harmonise_action)
  dat <- tryCatch(harmonise(exp_dat), error = function(err) {
    # Fall back to harmonising one exposure at a time
    map_dfr(split(exp_dat, exp_dat$exposure), function(x) tryCatch(harmonise(x), error = function(e) {
      message("skipped ", x$exposure[1], " / ", outcome_name, ": ", conditionMessage(e))
      NULL
    }))
  })
  if (nrow(dat) == 0) return(NULL)
  dat <- filter(dat, mr_keep)
  if (nrow(dat) == 0) return(NULL)

  snps <- dat %>% group_by(exposure) %>% summarise(snps = paste(SNP, collapse = ";"))
  TwoSampleMR::mr(dat, method_list = cfg$mr$methods) %>%
    as_tibble() %>%
    select(exposure, outcome, method, nsnp, b, se, pval) %>%
    left_join(snps, by = "exposure")
}
