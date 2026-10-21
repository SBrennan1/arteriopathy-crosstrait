# Sensitivity of the ASSET and CPASSOC loci to inflation: Z-scores divided by the square
# root of each trait's LDSC intercept (if > 1), both methods rerun within +/- 500 kb of
# each locus, and each locus rechecked for significance and cross-trait support.

source("R/setup.R")
library(ASSET)
source(cfg$cpassoc_functions)
select <- dplyr::select
out <- res_dir("sensitivity")
w <- 500e3

h2 <- read_tsv(file.path(cfg$results_dir, "ldsc/h2_liability.tsv"), show_col_types = FALSE)
z_sd <- set_names(sqrt(pmax(1, set_names(h2$intercept, h2$trait)[traits])), traits)
print(round(z_sd^2, 3))

zm <- load_z_matrix(ref_trait = "CeAD")
pos <- snp_positions(rownames(zm$z))

asset_loci <- read_tsv(file.path(cfg$results_dir, "asset/asset_loci.tsv"), show_col_types = FALSE)
cp_loci <- read_tsv(file.path(cfg$results_dir, "cpassoc/cpassoc_loci.tsv"), show_col_types = FALSE)
win <- bind_rows(asset_loci %>% select(chr, pos), cp_loci %>% select(chr, pos)) %>%
  transmute(chr = as.integer(chr), start = pos - w, end = pos + w)
keep <- rows_in_windows(pos %>% filter(!is.na(chr), !is.na(pos)), win)$rsid
keep <- intersect(rownames(zm$z), keep)
message(length(keep), " SNPs in the locus windows")

z <- sweep(zm$z[keep, ], 2, z_sd, "/")
p <- 2 * pnorm(-abs(z))
corr <- read_gcov_int()

# ASSET as in 01_asset.R
fit <- h.traits(snp.vars = keep, traits.lab = traits, beta.hat = z,
                sigma.hat = matrix(1, length(keep), length(traits), dimnames = list(NULL, traits)),
                cor = corr, ncase = map_dbl(cfg$traits, "ncase"), ncntl = map_dbl(cfg$traits, "ncontrol"),
                cor.numr = FALSE, side = 2, search = 2)
asset_p <- set_names(fit$Subset.2sided$pval, fit$snp.vars)

# CPASSOC SHet as in 02_cpassoc.R
para <- EstimateGamma(N = 1e4, SampleSize = zm$n, CorrMatrix = corr, correct = 1,
                      startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = TRUE)
s_het <- SHet(X = z, SampleSize = zm$n, CorrMatrix = corr)
shet_p <- set_names(pgamma(s_het - para[3], shape = para[1], scale = para[2], lower.tail = FALSE), keep)

# Best SNP in each locus window after rescaling, with cross-trait support at that SNP
# or at the unscaled lead, whichever is broader
support <- function(snp) {
  if (length(snp) == 0 || !snp %in% keep) return(nominal_support(matrix(NA, 1, length(traits), dimnames = list(NULL, traits)),
                                                                 matrix(NA, 1, length(traits), dimnames = list(NULL, traits))))
  nominal_support(z[snp, , drop = FALSE], p[snp, , drop = FALSE])
}
traits_of <- function(s) paste(na.omit(c(s$traits_pos, s$traits_neg)), collapse = ",")

recheck <- function(loci, stat, p_col) {
  pos_k <- pos %>% filter(rsid %in% keep)
  pmap_dfr(loci %>% select(rsid, chr, pos, traits_pos, traits_neg, all_of(p_col)),
           function(rsid, chr, pos, traits_pos, traits_neg, ...) {
             s <- pos_k$rsid[pos_k$chr == chr & abs(pos_k$pos - pos) <= w]
             best <- if (length(s)) s[which.min(stat[s])] else character()
             at_best <- support(best)
             at_lead <- support(rsid)
             use_lead <- coalesce(at_lead$n_nominal, 0) > coalesce(at_best$n_nominal, 0)
             tibble(locus_rsid = rsid, chr, pos, p_unscaled = list(...)[[1]],
                    traits_unscaled = paste(na.omit(c(traits_pos, traits_neg)), collapse = ","),
                    best_snp = if (length(best)) best else NA_character_,
                    p_rescaled = if (length(best)) unname(stat[best]) else NA_real_,
                    traits_rescaled = traits_of(if (use_lead) at_lead else at_best),
                    support_at = if (use_lead) "unscaled lead" else "best SNP",
                    n_nominal = max(coalesce(at_lead$n_nominal, 0), coalesce(at_best$n_nominal, 0)))
           }) %>%
    mutate(gws = p_rescaled < cfg$clump$p,
           cross_trait = n_nominal >= 2,
           retained = gws & cross_trait,
           aa_unscaled = str_detect(traits_unscaled, "\\bAA\\b"),
           aa_rescaled = str_detect(traits_rescaled, "\\bAA\\b"))
}

res_asset <- recheck(asset_loci, asset_p, "asset_p")
res_cp <- recheck(cp_loci, shet_p, "P_SHet")
write_tsv(res_asset, file.path(out, "intercept_asset.tsv"))
write_tsv(res_cp, file.path(out, "intercept_cpassoc.tsv"))

summ <- function(d) c(loci = nrow(d), retained = sum(d$retained, na.rm = TRUE), gws = sum(d$gws, na.rm = TRUE),
                      aa_loci = sum(d$aa_unscaled, na.rm = TRUE),
                      aa_still_aa = sum(d$aa_unscaled & d$aa_rescaled & d$retained, na.rm = TRUE))
print(rbind(ASSET = summ(res_asset), CPASSOC = summ(res_cp)))
write_counts("intercept_sensitivity", c(set_names(summ(res_asset), paste0("asset_", names(summ(res_asset)))),
                                        set_names(summ(res_cp), paste0("cpassoc_", names(summ(res_cp))))))
