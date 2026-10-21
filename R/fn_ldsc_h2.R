# SNP heritability with ldsc --h2 (liability scale when prevalences are given)

ldsc_h2 <- function(munged, out_prefix, samp_prev = NULL, pop_prev = NULL) {
  log <- paste0(out_prefix, ".log")
  # Rerun if the log used other prevalences
  if (file.exists(log) && !is.null(samp_prev)) {
    x <- readLines(log)
    logged <- function(flag) suppressWarnings(as.numeric(str_match(x, paste0("^", flag, " (\\S+)"))[, 2]))
    same <- function(flag, v) { l <- logged(flag); l <- l[!is.na(l)]; length(l) == 1 && isTRUE(all.equal(l, as.numeric(v))) }
    if (!same("--samp-prev", samp_prev) || !same("--pop-prev", pop_prev)) file.remove(log)
  }
  if (!file.exists(log)) {
    run_ldsc(c("--h2", munged, "--ref-ld-chr", cfg$ldsc$ref_ld, "--w-ld-chr", cfg$ldsc$w_ld,
               if (!is.null(samp_prev)) c("--samp-prev", samp_prev, "--pop-prev", pop_prev),
               "--out", out_prefix), log)
  }
  parse_h2_log(log)
}

# Estimate and standard error from a log line such as "Intercept: 1.0214 (0.0063)"
ldsc_value <- function(x, label) {
  num <- "NA|-?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
  m <- str_match(x, paste0("^", label, "\\s*(", num, ")(?:\\s*\\((", num, ")\\))?"))
  m <- m[!is.na(m[, 1]), , drop = FALSE]
  if (nrow(m) == 0) return(c(NA_real_, NA_real_))
  suppressWarnings(as.numeric(m[1, 2:3]))
}

parse_h2_log <- function(log) {
  x <- readLines(log)
  obs <- ldsc_value(x, "Total Observed scale h2:")
  liab <- ldsc_value(x, "Total Liability scale h2:")
  int <- ldsc_value(x, "Intercept:")
  rat <- ldsc_value(x, "Ratio:")
  # ldsc prints "Ratio < 0 (usually indicates GC correction)." instead of a value
  gc_corrected <- any(grepl("^Ratio < 0", x))
  n_snps <- as.numeric(str_match(x, "After merging with regression SNP LD, (\\d+) SNPs remain")[, 2])
  tibble(h2_obs = obs[1], h2_obs_se = obs[2], h2_liab = liab[1], h2_liab_se = liab[2],
         n_snps = n_snps[!is.na(n_snps)][1],
         lambda_gc = ldsc_value(x, "Lambda GC:")[1], mean_chi2 = ldsc_value(x, "Mean Chi\\^2:")[1],
         intercept = int[1], intercept_se = int[2], ratio = rat[1], ratio_se = rat[2],
         note = if (gc_corrected) "GC-corrected (Ratio < 0)" else NA_character_)
}
