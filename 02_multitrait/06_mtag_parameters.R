# Parameters of each pairwise MTAG run (N, mean chi2, GWAS-equivalent N, maxFDR and
# genetic correlation), parsed from the MTAG logs.

source("R/setup.R")
out <- res_dir("mtag")

num <- "-?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?"

parse_mtag_log <- function(log) {
  x <- readLines(log, encoding = "UTF-8")

  max_fdr <- map_dbl(1:2, function(i) {
    line <- grep(paste0("Max FDR of Trait ", i, ":"), x, value = TRUE)
    as.numeric(str_match(line[1], paste0("Trait ", i, ": (", num, ")"))[, 2])
  })

  # Summary rows: last six fields are SNPs, N (max), N (mean), GWAS chi2, MTAG chi2,
  # GWAS-equivalent N
  i <- grep("Summary of MTAG results:", x)[1]
  rows <- x[(i + 1):(i + 10)]
  rows <- rows[grepl("^\\d", rows)][1:2]
  summ <- map(rows, ~ tail(as.numeric(str_extract_all(.x, num)[[1]]), 6))

  # Genetic correlation from the first Omega correlation block
  j <- grep("^\\(Correlation\\):", x)[1]
  rg <- if (is.na(j)) NA_real_ else as.numeric(str_extract_all(x[j + 1], num)[[1]])[2]

  tibble(trait_index = 1:2,
         n = map_dbl(summ, 2),
         gwas_chi2 = map_dbl(summ, 4),
         mtag_chi2 = map_dbl(summ, 5),
         gwas_equiv_n = map_dbl(summ, 6),
         max_fdr = max_fdr,
         rg_mtag = rg)
}

params <- imap_dfr(cfg$mtag_logs, function(log, pair) {
  parse_mtag_log(file.path(cfg$mtag$dir, log)) %>%
    mutate(mtag_pair = pair, trait = str_split_1(pair, "-")[trait_index])
}) %>%
  mutate(gwas_chi2 = round(gwas_chi2, 3),
         mtag_chi2 = round(mtag_chi2, 3),
         max_fdr = round(max_fdr, 3),
         rg_mtag = round(rg_mtag, 3),
         chi2_change = round(mtag_chi2 - gwas_chi2, 3)) %>%
  select(mtag_pair, trait, n, rg_mtag, gwas_chi2, mtag_chi2, chi2_change, gwas_equiv_n, max_fdr)

write_tsv(params, file.path(out, "mtag_parameters.tsv"))

write_counts("mtag_parameters",
             setNames(params$max_fdr, paste0("max_fdr_", params$mtag_pair, "_", params$trait)))
