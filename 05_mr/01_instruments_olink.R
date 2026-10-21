# UKB-PPP (Olink) cis-pQTL instruments: P < 1.7e-11, F >= 10, clumped per protein
# (r2 < 0.001 within 250 kb, 1000G EUR).

source("R/setup.R")
out <- res_dir("mr")
ol <- cfg$mr$olink

olink <- vroom(ol$raw, show_col_types = FALSE)
stopifnot(all(c("variant", "Log10_pval_gc_cor_unadj", "cis/trans", "beta", "se", "eaf",
                "rsid", "UniProt", "protein") %in% names(olink)))

# variant is CHR:POS:A0:A1:imp:v1 with GRCh37 positions; beta and eaf refer to A1
v <- str_split_fixed(olink$variant, ":", 5)

cis <- olink %>%
  mutate(chr = v[, 1], pos = as.integer(v[, 2]), other_allele = v[, 3], effect_allele = v[, 4],
         across(c(eaf, beta, se, Log10_pval_gc_cor_unadj), as.numeric),
         eaf = if_else(eaf < 0 | eaf > 1, NA_real_, eaf)) %>%
  filter(tolower(`cis/trans`) == "cis") %>%
  mutate(
    lp = Log10_pval_gc_cor_unadj,
    pval = case_when(is.finite(lp) & lp >= 0 ~ 10^-lp,
                     is.infinite(lp) & lp > 0 ~ 0),
    pval = if_else(pval == 0, .Machine$double.xmin, pval),
    z = if_else(is.finite(se) & se != 0, beta / se, NA_real_),
    z = if_else(is.finite(z), z, NA_real_),
    F_statistic = z^2
  )

strong <- cis %>%
  filter(!is.na(pval), pval < ol$p, !is.na(F_statistic), F_statistic >= ol$f_min)

# Clumping ----
clump_protein <- function(d) {
  d <- filter(d, !is.na(rsid), !is.na(pval))
  if (nrow(d) <= 1) return(d)
  kept <- tryCatch(clump_local(d, p = 1, r2 = ol$clump_r2, kb = ol$clump_kb),
                   error = function(e) NULL)
  semi_join(d, kept %||% tibble(rsid = character()), by = "rsid")
}

instruments <- split(strong, strong$UniProt) %>%
  map_dfr(clump_protein) %>%
  transmute(protein, UniProt, rsid, chr, pos, effect_allele, other_allele, eaf, beta, se,
            pval, F_statistic, variant) %>%
  add_hg38()

write_tsv(instruments, file.path(out, "instruments_olink.tsv"))

write_counts("mr_instruments", c(
  olink_cis_variants = nrow(cis),
  olink_strong_variants = nrow(strong),
  olink_instruments = nrow(instruments),
  olink_proteins = n_distinct(instruments$protein)
))
