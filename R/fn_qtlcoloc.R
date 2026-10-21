# Helpers for the molecular QTL colocalisation stage (06_qtl_coloc)

# Runs an expression without coloc/TwoSampleMR console output
run_quiet <- function(x) {
  utils::capture.output(x <- suppressMessages(x))
  x
}

# +/- window around each protein's gene body (GRCh37), from cfg$qtl_coloc$regions or
# the Ensembl GTF
protein_regions <- function(proteins, kb = cfg$qtl_coloc$window_kb) {
  fixed <- read_tsv(cfg$qtl_coloc$regions, show_col_types = FALSE,
                    col_types = cols(protein = "c", chr = "c", start = "i", end = "i"))
  todo <- setdiff(proteins, fixed$protein)
  gtf <- if (length(todo)) {
    # any biotype: some measured proteins (e.g. SPON1) are not protein_coding in Ensembl 87
    load_genes(protein_coding = FALSE) %>%
      filter(gene_name %in% todo, chr %in% as.character(1:22)) %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      transmute(protein = gene_name, chr, start = as.integer(start), end = as.integer(end))
  }
  tibble(protein = unique(proteins)) %>%
    left_join(bind_rows(fixed, gtf), by = "protein") %>%
    transmute(protein, chr = as.integer(chr),
              start = pmax(1L, start - kb * 1000L), end = end + kb * 1000L)
}

# Rows of d (chr, pos) falling in any of the windows (chr, start, end)
rows_in_windows <- function(d, win) {
  hit <- as.data.table(d)[as.data.table(win), on = .(chr, pos >= start, pos <= end),
                          nomatch = 0, which = TRUE]
  d[sort(unique(hit)), ]
}

# pQTL exposures ----

# UKB-PPP rsID/GRCh37 map, one file per chromosome, read once per session
olink_rsid_map <- local({
  cache <- list()
  function(chr) {
    key <- as.character(chr)
    if (is.null(cache[[key]])) {
      f <- list.files(cfg$qtl_coloc$olink$rsid_map_dir, pattern = paste0("_chr", chr, "_"),
                      full.names = TRUE)
      if (length(f) == 0) stop("no UKB-PPP rsID map for chr", chr)
      cache[[key]] <<- vroom(f, show_col_types = FALSE, col_select = c(ID, rsid, POS19))
    }
    cache[[key]]
  }
})

# UKB-PPP discovery summary statistics (GRCh38) for the cis chromosome, put on
# GRCh37 rsIDs/positions through the UKB-PPP map
read_olink_region <- function(protein, region) {
  dirs <- list.dirs(cfg$qtl_coloc$olink$dir, recursive = FALSE)
  dirs <- dirs[startsWith(toupper(basename(dirs)), paste0(toupper(protein), "_"))]
  f <- list.files(dirs, pattern = paste0("^discovery_chr", region$chr, "_"),
                  full.names = TRUE, recursive = TRUE)
  if (length(f) == 0) stop("no UKB-PPP file")
  # Proteins on several Olink panels: largest file
  f <- f[which.max(file.size(f))]

  vroom(f, show_col_types = FALSE,
        col_select = c(CHROM, ID, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P)) %>%
    filter(CHROM == region$chr) %>%
    inner_join(olink_rsid_map(region$chr), by = "ID") %>%
    transmute(SNP = rsid, chr = as.integer(CHROM), pos = as.integer(POS19),
              effect_allele = ALLELE1, other_allele = ALLELE0, eaf = A1FREQ,
              beta = BETA, se = SE, pval = 10^(-LOG10P)) %>%
    filter(pos >= region$start, pos <= region$end) %>%
    drop_na(SNP, pos, beta, se, effect_allele, other_allele) %>%
    mutate(pval = pmax(pval, 1e-300)) %>%
    distinct(SNP, .keep_all = TRUE)
}

# deCODE regional files ({PROTEIN}.csv, already on GRCh37)
read_decode_region <- function(protein, region) {
  files <- list.files(cfg$qtl_coloc$decode$dir, pattern = "\\.csv$", full.names = TRUE)
  f <- files[toupper(tools::file_path_sans_ext(basename(files))) == toupper(protein)]
  if (length(f) == 0) stop("no deCODE file")

  vroom(f[1], show_col_types = FALSE) %>%
    transmute(SNP = rsids, chr = as.integer(str_remove(as.character(chr), "^chr")),
              pos = as.integer(pos), effect_allele = toupper(effectAllele),
              other_allele = toupper(otherAllele), eaf = as.numeric(ImpMAF),
              beta = as.numeric(Beta), se = as.numeric(SE), pval = as.numeric(Pval)) %>%
    separate_rows(SNP, sep = "[,;\\s]+") %>%
    mutate(SNP = trimws(SNP)) %>%
    filter(str_detect(SNP, "^rs\\d+$"), chr == region$chr,
           pos >= region$start, pos <= region$end,
           !is.na(beta), !is.na(se), !is.na(pval)) %>%
    mutate(pval = pmax(pval, 1e-300)) %>%
    distinct(SNP, .keep_all = TRUE)
}

# pQTL coloc ----

harmonise_qtl <- function(exposure, outcome) {
  fmt <- function(d, type) {
    run_quiet(TwoSampleMR::format_data(
      d, type = type, snp_col = "SNP", chr_col = "chr", pos_col = "pos",
      beta_col = "beta", se_col = "se", pval_col = "pval", eaf_col = "eaf",
      effect_allele_col = "effect_allele", other_allele_col = "other_allele"))
  }
  exp <- fmt(exposure, "exposure")
  out <- fmt(rename(outcome, SNP = rsid), "outcome")
  run_quiet(TwoSampleMR::harmonise_data(exp, out, action = 2)) %>%
    filter(!is.na(beta.exposure), !is.na(beta.outcome)) %>%
    distinct(SNP, .keep_all = TRUE)
}

# Signed LD (r) from the per-chromosome 1000G EUR panel
ld_r_chr <- function(rsids, chr) {
  R <- ieugwasr::ld_matrix(rsids, with_alleles = FALSE, bfile = file.path(cfg$ld_ref_chr, chr),
                           plink_bin = plink_bin())
  keep <- !duplicated(rownames(R))
  R[keep, keep, drop = FALSE]
}

# Drops the later SNP of any pair with |r| >= thr
pqtl_prune_ld <- function(R, thr) {
  keep <- rep(TRUE, nrow(R))
  for (i in seq_len(nrow(R))) {
    if (!keep[i]) next
    j <- which(abs(R[i, ]) >= thr)
    keep[j[j > i]] <- FALSE
  }
  keep
}

# Number of SuSiE effects: 1 + (number of |z| > 5) / 60, capped at 10
pqtl_susie_L <- function(z) as.integer(max(1L, min(10L, round(sum(abs(z) > 5) / 60) + 1L)))

# Credible sets whose lead SNPs are in high LD are merged (keeping the set with the
# larger maximum alpha), so one signal is not split across several sets
pqtl_merge_cs <- function(fit, R, r2_thr) {
  idx <- fit$sets$cs_index
  if (length(idx) < 2) return(fit)
  lead <- vapply(idx, function(k) colnames(R)[which.max(fit$alpha[k, ])], character(1))
  pip <- vapply(idx, function(k) max(fit$alpha[k, ], na.rm = TRUE), numeric(1))
  keep <- rep(TRUE, length(idx))
  for (i in seq_along(idx)[-length(idx)]) {
    if (!keep[i]) next
    for (j in (i + 1):length(idx)) {
      if (keep[j] && R[lead[i], lead[j]]^2 >= r2_thr) {
        if (pip[i] >= pip[j]) keep[j] <- FALSE else keep[i] <- FALSE
      }
    }
  }
  fit$sets$cs <- fit$sets$cs[keep]
  fit$sets$cs_index <- idx[keep]
  fit$sets$coverage <- fit$sets$coverage[keep]
  if (!is.null(fit$sets$purity)) fit$sets$purity <- fit$sets$purity[keep, , drop = FALSE]
  fit
}

# SuSiE on both traits with 1000G LD, then coloc.susie
pqtl_susie_coloc <- function(dat, chr, n_exp, n_out, p = cfg$qtl_coloc) {
  R <- ld_r_chr(dat$SNP, chr)
  common <- intersect(dat$SNP, rownames(R))
  if (length(common) < p$min_snps) stop("too few SNPs in the LD panel")
  R <- R[common, common, drop = FALSE]
  if (any(!is.finite(R))) stop("LD matrix has undefined correlations")
  dat <- dat[match(common, dat$SNP), ]

  keep <- pqtl_prune_ld(R, p$ld_prune_r)
  R <- R[keep, keep, drop = FALSE]
  dat <- dat[keep, ]
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R <- (1 - p$ridge) * R + p$ridge * diag(nrow(R))
  if (inherits(try(chol(R), silent = TRUE), "try-error")) {
    R <- as.matrix(Matrix::nearPD(R, keepDiag = TRUE)$mat)
  }

  z_exp <- dat$beta.exposure / dat$se.exposure
  z_out <- dat$beta.outcome / dat$se.outcome
  good <- is.finite(z_exp) & is.finite(z_out)
  if (sum(good) < p$min_snps) stop("too few finite z-scores")
  dat <- dat[good, ]
  R <- R[good, good, drop = FALSE]

  # z-score SuSiE with residual variance fixed at 1 (inverse-normalised proteins)
  fit <- function(z, n) {
    f <- susieR::susie_rss(z = z, R = R, n = n, L = pqtl_susie_L(z),
                           estimate_residual_variance = FALSE, residual_variance = 1,
                           max_iter = p$susie_max_iter)
    pqtl_merge_cs(f, R, p$cs_merge_r2)
  }
  fit_exp <- fit(z_exp[good], n_exp)
  fit_out <- fit(z_out[good], n_out)

  s <- tryCatch(
    run_quiet(coloc::coloc.susie(fit_exp, fit_out, p1 = p$p1, p2 = p$p2, p12 = p$p12))$summary,
    error = function(e) NULL)
  ok <- is.data.frame(s) && nrow(s) > 0
  list(dat = dat,
       susie_pph3 = if (ok) max(s$PP.H3.abf) else NA_real_,
       susie_pph4 = if (ok) max(s$PP.H4.abf) else NA_real_,
       susie_hit1 = if (ok) as.character(s$hit1[which.max(s$PP.H4.abf)]) else NA_character_)
}

# Protein (quantitative, sdY = 1) vs case-control outcome
abf_pqtl <- function(dat, n_exp, ncase, ncontrol, p = cfg$qtl_coloc) {
  n_out <- ncase + ncontrol
  run_quiet(coloc::coloc.abf(
    list(snp = dat$SNP, beta = dat$beta.exposure, varbeta = dat$se.exposure^2,
         type = "quant", sdY = 1, N = n_exp),
    list(snp = dat$SNP, beta = dat$beta.outcome, varbeta = dat$se.outcome^2,
         type = "cc", N = n_out, s = ncase / n_out),
    p1 = p$p1, p2 = p$p2, p12 = p$p12))$summary
}

coloc_pqtl <- function(exposure, outcome, region, n_exp, ncase, ncontrol, p = cfg$qtl_coloc) {
  if (nrow(exposure) < p$min_snps) stop("too few exposure SNPs")
  if (nrow(outcome) < p$min_snps) stop("too few outcome SNPs")
  dat <- harmonise_qtl(exposure, outcome)
  if (nrow(dat) < p$min_snps) stop("too few harmonised SNPs")

  n_eff <- 4 / (1 / ncase + 1 / ncontrol)
  ss <- tryCatch(pqtl_susie_coloc(dat, region$chr, n_exp, n_eff, p = p), error = function(e) NULL)
  if (!is.null(ss)) dat <- ss$dat
  abf <- abf_pqtl(dat, n_exp, ncase, ncontrol, p)
  lead_exp <- dat$SNP[which.max(abs(dat$beta.exposure / dat$se.exposure))]

  tibble(chr = region$chr, start = region$start, end = region$end, nsnps = nrow(dat),
         susie_pph3 = ss$susie_pph3 %||% NA_real_,
         susie_pph4 = ss$susie_pph4 %||% NA_real_,
         abf_pph3 = unname(abf["PP.H3.abf"]),
         abf_pph4 = unname(abf["PP.H4.abf"]),
         lead_exp = lead_exp,
         lead_out = dat$SNP[which.max(abs(dat$beta.outcome / dat$se.outcome))],
         coloc_snp = coalesce(ss$susie_hit1 %||% NA_character_, lead_exp),
         method = if (is.null(ss)) "ABF only" else "SuSiE+ABF")
}

# eQTL coloc ----

lift_eqtl_hg19 <- function(chr, pos, chain) {
  out <- rep(NA_integer_, length(pos))
  ok <- which(!is.na(chr) & !is.na(pos))
  gr <- GenomicRanges::GRanges(paste0("chr", chr[ok]),
                               IRanges::IRanges(as.integer(pos[ok]), width = 1))
  lo <- rtracklayer::liftOver(gr, chain)
  one <- which(lengths(lo) == 1)
  hit <- unlist(lo[one])
  same <- as.character(GenomicRanges::seqnames(hit)) == paste0("chr", chr[ok][one])
  out[ok[one[same]]] <- GenomicRanges::start(hit)[same]
  out
}

# GTEx v10 regional subset for one tissue-gene pair, lifted to GRCh37
read_eqtl_subset <- function(tissue, gene_id, files, chain, p = cfg$qtl_coloc$eqtl) {
  tissue_file <- p$tissue_files[[tissue]] %||% tissue
  f <- grep(paste0(tissue_file, "\\.", sub("\\.\\d+$", "", gene_id), "\\."), files, value = TRUE)
  if (length(f) == 0) stop("no eQTL subset file")
  d <- fread(f[1]) %>%
    as_tibble() %>%
    mutate(chr = as.integer(str_remove(as.character(chr), "^chr")),
           across(c(pos, af, pval_nominal, slope, slope_se, ma_samples), as.numeric),
           pos37 = lift_eqtl_hg19(chr, pos, chain)) %>%
    filter(!is.na(pos37))
  if (nrow(d) < p$min_snps) stop("too few variants lifted to GRCh37")
  d
}

# eQTL (quantitative, sdY estimated from MAF and N) vs case-control outcome, matched
# by GRCh37 position. gwas is a data.table with chr, pos, alleles, beta, se.
coloc_eqtl <- function(eq, gwas, ncase, ncontrol, p = cfg$qtl_coloc$eqtl) {
  need <- function(d, what) {
    if (nrow(d) < p$min_snps) stop("too few SNPs ", what)
    d
  }
  # eQTL N: median ma_samples
  n_eqtl <- median(eq$ma_samples, na.rm = TRUE)

  # GTEx variant_id is chr_pos_ref_alt_b38 and slope is the effect of alt
  alleles <- str_split_fixed(eq$variant_id, "_", 5)
  eq <- eq %>%
    transmute(chr, pos = as.integer(pos37), ref = alleles[, 3], alt = alleles[, 4],
              af, pval_nominal, slope, slope_se) %>%
    as.data.table()

  m <- gwas[eq, on = c("chr", "pos"), nomatch = NULL] %>%
    as_tibble() %>%
    need("overlapping the GWAS") %>%
    mutate(fwd = alt == effect_allele & ref == other_allele,
           rev = alt == other_allele & ref == effect_allele,
           slope = if_else(rev, -slope, slope),
           af = if_else(rev, 1 - af, af)) %>%
    filter(fwd | rev) %>%
    need("after allele matching") %>%
    filter(!is.na(slope), !is.na(slope_se), !is.na(beta), !is.na(se)) %>%
    need("with complete statistics") %>%
    # several GRCh38 variants can land on one GRCh37 position: keep the strongest eQTL
    arrange(pval_nominal) %>%
    distinct(chr, pos, .keep_all = TRUE) %>%
    need("after removing duplicate positions")

  snp <- paste(m$chr, m$pos, sep = ":")
  r <- run_quiet(coloc::coloc.abf(
    list(snp = snp, beta = m$slope, varbeta = m$slope_se^2, MAF = m$af, N = n_eqtl, type = "quant"),
    list(snp = snp, beta = m$beta, varbeta = m$se^2, N = ncase + ncontrol,
         s = ncase / (ncase + ncontrol), type = "cc"),
    p1 = p$p1, p2 = p$p2, p12 = p$p12))$summary
  tibble(nsnps = unname(r["nsnps"]),
         PP.H0 = unname(r["PP.H0.abf"]), PP.H1 = unname(r["PP.H1.abf"]),
         PP.H2 = unname(r["PP.H2.abf"]), PP.H3 = unname(r["PP.H3.abf"]),
         PP.H4 = unname(r["PP.H4.abf"]), method = "ABF")
}

# Colocalization evidence class (PP.H4 rounded to 3 dp)
qtl_coloc_evidence <- function(pph4, strong = cfg$qtl_coloc$strong,
                               suggestive = cfg$qtl_coloc$suggestive) {
  x <- round(pph4, 3)
  case_when(is.na(x) ~ "Not tested",
            x >= strong ~ "Strong",
            x >= suggestive ~ "Suggestive",
            TRUE ~ "None")
}
