run_ldsc <- function(args, log) {
  out <- system2(cfg$ldsc$ldsc, args, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status")) || !file.exists(log)) {
    stop("ldsc failed:\n", paste(tail(out, 20), collapse = "\n"))
  }
  invisible(log)
}

# hm3 = TRUE: LDSC input (HapMap3 SNPs, default filters).
# hm3 = FALSE: genome-wide input for ASSET, CPASSOC and MTAG, without the MAF filter.
munge <- function(infile, out_prefix, hm3 = TRUE) {
  if (file.exists(paste0(out_prefix, ".sumstats.gz"))) return(invisible())
  # Use eaf only where available (CeAD has none)
  use_frq <- hm3 &&
    any(!is.na(vroom(infile, n_max = 1e5, col_select = "eaf", show_col_types = FALSE)$eaf))
  # munge_sumstats.py (Python 3) reads an uncompressed copy
  if (grepl("\\.gz$", infile)) {
    plain <- tempfile(fileext = ".tsv")
    on.exit(unlink(plain), add = TRUE)
    if (system2("gzip", c("-dc", shQuote(infile)), stdout = plain) != 0) stop("could not decompress ", infile)
    infile <- plain
  }
  out <- system2(cfg$ldsc$munge, c(
    "--sumstats", infile, "--out", out_prefix,
    "--snp", "rsid", "--a1", "effect_allele", "--a2", "other_allele",
    "--p", "pval", "--signed-sumstats", "beta,0", "--N-col", "N",
    if (use_frq) c("--frq", "eaf") else c("--ignore", "eaf"),
    if (hm3) c("--merge-alleles", cfg$ldsc$hm3),
    "--chunksize", "500000"
  ), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status"))) {
    stop("munge failed for ", infile, ":\n", paste(tail(out, 15), collapse = "\n"))
  }
}

ldsc_rg <- function(files, out_prefix) {
  log <- paste0(out_prefix, ".log")
  # Reuse an existing log only if it covers the same files
  if (file.exists(log)) {
    old <- tryCatch(parse_rg_log(log), error = function(e) NULL)
    want <- sub("\\.sumstats\\.gz$", "", basename(files[-1]))
    if (is.null(old) || !setequal(old$p2, want)) file.remove(log)
  }
  if (!file.exists(log)) {
    run_ldsc(c("--rg", paste(files, collapse = ","),
               "--ref-ld-chr", cfg$ldsc$ref_ld, "--w-ld-chr", cfg$ldsc$w_ld,
               "--out", out_prefix), log)
  }
  parse_rg_log(log)
}

# Reads the summary table at the end of an ldsc --rg log
parse_rg_log <- function(log) {
  x <- readLines(log)
  i <- grep("^Summary of Genetic Correlation Results", x)
  if (length(i) == 0) stop("no rg summary in ", log)
  tab <- x[(i + 1):length(x)]
  tab <- tab[nzchar(trimws(tab)) & !grepl("^Analysis finished|^Total time", tab)]
  read.table(text = tab, header = TRUE, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    mutate(across(-c(p1, p2), as.numeric),
           p1 = sub("\\.sumstats\\.gz$", "", basename(p1)),
           p2 = sub("\\.sumstats\\.gz$", "", basename(p2)))
}
