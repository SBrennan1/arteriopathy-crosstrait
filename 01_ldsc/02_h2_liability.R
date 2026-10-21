# Observed- and liability-scale SNP heritability (LDSC), using the prevalences in
# config (ldsc_h2).

source("R/setup.R")
out <- res_dir("ldsc", "h2")

# Run with and without prevalences to obtain both scales
h2 <- map_dfr(traits, function(t) {
  prev <- cfg$ldsc_h2$prevalence[[t]]
  f <- file.path(cfg$munged_dir, paste0(t, ".sumstats.gz"))
  liab <- ldsc_h2(f, file.path(out, paste0(t, "_h2_liab")), samp_prev = prev$samp, pop_prev = prev$pop)
  ldsc_h2(f, file.path(out, paste0(t, "_h2"))) %>%
    mutate(h2_liab = liab$h2_liab, h2_liab_se = liab$h2_liab_se) %>%
    mutate(trait = t, samp_prev = prev$samp, pop_prev = prev$pop, .before = 1)
})

write_tsv(h2, file.path(res_dir("ldsc"), "h2_liability.tsv"))

write_counts("ldsc_h2", c(
  setNames(round(h2$h2_liab, 4), paste0("h2_liab_", h2$trait)),
  setNames(round(h2$h2_liab_se, 4), paste0("h2_liab_se_", h2$trait))
))
