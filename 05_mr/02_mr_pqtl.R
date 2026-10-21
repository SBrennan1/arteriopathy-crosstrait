# Proteome-wide MR of cis-pQTL instruments (UKB-PPP Olink, deCODE SomaScan) on the
# five arteriopathies.

source("R/setup.R")
out <- res_dir("mr")

read_pqtl <- function(path) {
  d <- vroom(path, show_col_types = FALSE,
             col_types = cols(rsid = "c", effect_allele = "c", other_allele = "c"))
  name_col <- intersect(c("exposure_name", "protein", "exposure"), names(d))[1]
  d <- d %>%
    transmute(exposure = .data[[name_col]], rsid = as.character(rsid), effect_allele,
              other_allele, beta = as.numeric(beta), se = as.numeric(se),
              eaf = suppressWarnings(as.numeric(eaf)), pval = as.numeric(pval))
  if (any(is.na(d$se))) stop("missing standard errors in ", path)
  d
}

platforms <- names(cfg$mr$pqtl_instruments)
ivs <- map(set_names(platforms), ~ read_pqtl(cfg$mr$pqtl_instruments[[.x]]))
outcomes <- map(set_names(traits), read_mr_outcome, rsids = unique(unlist(map(ivs, "rsid"))))

mr_pqtl <- function(platform) {
  imap_dfr(outcomes, ~ run_mr(ivs[[platform]], .x, .y))
}

res <- map_dfr(set_names(platforms), mr_pqtl, .id = "platform")
write_tsv(res, file.path(out, "mr_pqtl_unadjusted.tsv"))
