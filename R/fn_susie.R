# SuSiE fine-mapping against UK Biobank LD (Weissbrod et al. 2020). The LD
# matrices are scipy sparse .npz files in overlapping 3 Mb chunks
# (chr{chr}_{start}_{end}.npz, SNP table in the matching .gz), so reading them
# needs python with numpy and scipy, reached through reticulate.

ukbb_ld_init <- function(python = cfg$ukbb_ld$python) {
  if (!is.null(python)) reticulate::use_python(python, required = TRUE)
  if (!reticulate::py_module_available("scipy")) {
    stop("Reading the UKBB LD chunks needs python with numpy and scipy (ukbb_ld$python in config.yml)")
  }
  reticulate::py_run_string(paste(
    "import numpy as np",
    "import scipy.sparse as sp",
    "def load_and_subset_npz(npz_path, indices):",
    "    idx = np.array(indices)",
    "    return sp.load_npz(npz_path).tocsc()[:, idx].tocsr()[idx, :].toarray()",
    sep = "\n"
  ))
  invisible(TRUE)
}

ukbb_ld_catalogue <- function(dir = cfg$ukbb_ld$dir) {
  f <- list.files(dir, pattern = "^chr[0-9]+_[0-9]+_[0-9]+\\.npz$")
  m <- str_match(f, "^chr([0-9]+)_([0-9]+)_([0-9]+)\\.npz$")
  tibble(chr = as.integer(m[, 2]), start = as.numeric(m[, 3]), end = as.numeric(m[, 4]),
         npz = file.path(dir, f), snps = sub("\\.npz$", ".gz", file.path(dir, f)))
}

# LD among the panel SNPs in [start, end], taken from the single chunk that
# holds most of them. Returns list(R, snps) or NULL.
load_ukbb_ld <- function(catalogue, chr, start, end) {
  chunks <- catalogue[catalogue$chr == chr & catalogue$start <= end & catalogue$end >= start, ]
  best <- NULL
  best_n <- 1
  for (i in seq_len(nrow(chunks))) {
    snps <- fread(chunks$snps[i], header = TRUE)
    idx <- which(snps$position >= start & snps$position <= end)
    if (length(idx) <= best_n) next
    R <- reticulate::import_main(convert = TRUE)$load_and_subset_npz(chunks$npz[i], as.integer(idx - 1L))
    # chunks store one triangle
    if (!isSymmetric(R, tol = 1e-6)) {
      R <- R + t(R)
      diag(R) <- diag(R) / 2
    }
    diag(R) <- 1
    dimnames(R) <- list(snps$rsid[idx], snps$rsid[idx])
    best <- list(R = R, snps = snps[idx])
    best_n <- length(idx)
  }
  best
}

# Keep GWAS SNPs present in the LD panel with matching alleles; beta is put on
# the panel's allele1 so that z and R agree.
match_to_ld <- function(region, ld) {
  common <- intersect(region$rsid, ld$snps$rsid)
  if (length(common) < 2) return(NULL)
  region <- region[match(common, region$rsid), ]
  info <- ld$snps[match(common, ld$snps$rsid), ]
  a1 <- toupper(info$allele1)
  a2 <- toupper(info$allele2)
  same <- coalesce(region$effect_allele == a1 & region$other_allele == a2, FALSE)
  swap <- coalesce(region$effect_allele == a2 & region$other_allele == a1, FALSE) & !same
  region$beta[swap] <- -region$beta[swap]
  keep <- same | swap
  if (sum(keep) < 2) return(NULL)
  list(region = region[keep, ], R = ld$R[common, common, drop = FALSE][keep, keep, drop = FALSE],
       n_flipped = sum(swap))
}

# Shrink R towards the identity, trying each lambda in turn until it is positive
# definite (the last one is used regardless)
ridge_ld <- function(R, lambdas = cfg$susie$ridge) {
  diag(R) <- 1
  for (lambda in lambdas) {
    Rr <- (1 - lambda) * R + lambda * diag(nrow(R))
    diag(Rr) <- 1
    if (!inherits(try(chol(Rr), silent = TRUE), "try-error")) break
  }
  dimnames(Rr) <- dimnames(R)
  list(R = Rr, lambda = lambda)
}

# SuSiE (RSS) for one trait in one window. region: GWAS rows in the window
# (rsid, chr, pos, effect_allele, other_allele, beta, se, pval); n: GWAS sample size.
finemap_trait <- function(region, ld, n, p = cfg$susie) {
  region <- region %>%
    filter(!is.na(beta), !is.na(se), !is.na(effect_allele), !is.na(other_allele)) %>%
    distinct(rsid, .keep_all = TRUE)
  fail <- function(status, n_in_ld = NA_integer_) {
    list(status = status, n_tested = nrow(region), n_in_ld = n_in_ld, fit = NULL, has_cs = FALSE)
  }
  if (nrow(region) < p$min_snps) return(fail("too few GWAS SNPs"))

  m <- match_to_ld(region, ld)
  n_in_ld <- if (is.null(m)) 0L else nrow(m$region)
  if (n_in_ld < p$min_snps) return(fail("too few SNPs in LD panel", n_in_ld))

  r <- ridge_ld(m$R, p$ridge)
  z <- m$region$beta / m$region$se
  z[!is.finite(z)] <- 0
  fit <- tryCatch(
    susieR::susie_rss(z = z, R = r$R, n = n, L = p$L, coverage = p$coverage,
                      estimate_residual_variance = FALSE, residual_variance = 1,
                      max_iter = p$max_iter),
    error = function(e) NULL)
  if (is.null(fit)) return(fail("SuSiE error", n_in_ld))
  if (isFALSE(fit$converged)) return(fail("SuSiE did not converge", n_in_ld))

  cs <- fit$sets$cs
  idx <- unlist(cs)
  lead <- m$region %>% filter(is.finite(pval)) %>% slice_min(pval, n = 1, with_ties = FALSE)
  list(
    status = "ok", n_tested = nrow(region), n_in_ld = n_in_ld, n_flipped = m$n_flipped,
    ridge_lambda = r$lambda, fit = fit, has_cs = length(cs) > 0, n_cs = length(cs),
    cs_size = if (length(cs)) max(lengths(cs)) else NA_integer_,
    top_pip_snp = m$region$rsid[which.max(fit$pip)], top_pip = max(fit$pip),
    lead_snp = lead$rsid[1], lead_p = lead$pval[1],
    credible = if (length(cs)) {
      m$region[idx, ] %>%
        transmute(cs = rep(seq_along(cs), lengths(cs)), rsid, chr, pos,
                  effect_allele, other_allele, pip = unname(fit$pip[idx]))
    }
  )
}

finemap_summary <- function(f) {
  tibble(status = f$status, n_tested = f$n_tested, n_in_ld = f$n_in_ld,
         n_flipped = f$n_flipped %||% NA_integer_, ridge_lambda = f$ridge_lambda %||% NA_real_,
         n_cs = f$n_cs %||% NA_integer_, cs_size = f$cs_size %||% NA_integer_,
         top_pip_snp = f$top_pip_snp %||% NA_character_, top_pip = f$top_pip %||% NA_real_,
         lead_snp = f$lead_snp %||% NA_character_, lead_p = f$lead_p %||% NA_real_)
}

# coloc.susie between two fine-mapped traits. When both traits have several
# credible sets, the pair of sets with the highest PP.H4 is reported.
coloc_susie_pair <- function(f1, f2, p = cfg$susie) {
  co <- tryCatch(
    suppressMessages(coloc::coloc.susie(f1$fit, f2$fit, p1 = p$p1, p2 = p$p2, p12 = p$p12)),
    error = function(e) NULL)
  s <- co$summary
  if (is.null(s) || nrow(s) == 0 || !"PP.H4.abf" %in% names(s)) return(NULL)
  k <- nrow(s)
  as_tibble(s) %>%
    slice_max(PP.H4.abf, n = 1, with_ties = FALSE) %>%
    transmute(n_cs_pairs = k, hit1, hit2, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)
}

coloc_class <- function(h4, h3, p = cfg$susie) {
  case_when(
    h4 >= p$pph4_strong ~ "Strong",
    h4 >= p$pph4_moderate ~ "Moderate",
    h3 > 0.5 ~ "Distinct",
    TRUE ~ "None"
  )
}

# Traits with nominal support at an ASSET lead (positive direction first)
contributing_traits <- function(traits_pos, traits_neg) {
  map2(traits_pos, traits_neg, function(a, b) {
    x <- c(str_split_1(coalesce(a, ""), ","), str_split_1(coalesce(b, ""), ","))
    unique(x[x != ""])
  })
}
