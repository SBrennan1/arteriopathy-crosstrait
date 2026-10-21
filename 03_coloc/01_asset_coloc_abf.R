# Pairwise coloc.abf between the traits contributing to each ASSET pleiotropic locus
# (nominal support, P < 0.05), +/- 500 kb around the lead variant.

source("R/setup.R")
library(coloc)
out <- res_dir("coloc")

loci <- read_tsv(file.path(cfg$results_dir, "asset/asset_loci.tsv"), show_col_types = FALSE) %>%
  mutate(traits = map2(traits_pos, traits_neg, ~ na.omit(c(str_split_1(coalesce(.x, ""), ","),
                                                             str_split_1(coalesce(.y, ""), ","))))) %>%
  mutate(traits = map(traits, ~ .x[.x != ""]))

w <- cfg$coloc$window_kb * 1e3
regions <- loci %>% transmute(rsid, chr, start = pmax(1, pos - w), end = pos + w)

# Load each trait once, restricted to the locus windows
keep_regions <- function(d) {
  hit <- as.data.table(d)[as.data.table(regions), on = .(chr, pos >= start, pos <= end),
                          nomatch = 0, which = TRUE]
  d[sort(unique(hit)), ]
}

sumstats <- map(set_names(traits), function(t) {
  read_univar(t, cols = c("rsid", "chr", "pos", "beta", "se")) %>%
    filter(!is.na(beta), !is.na(se), se > 0) %>%
    keep_regions()
})

coloc_pair <- function(t1, t2, region) {
  pick_region <- function(t) sumstats[[t]] %>%
    filter(chr == region$chr, pos >= region$start, pos <= region$end) %>%
    distinct(rsid, .keep_all = TRUE)
  m <- inner_join(pick_region(t1), pick_region(t2), by = "rsid", suffix = c("1", "2")) %>%
    mutate(pos = pos1) %>%
    arrange(pos)

  res <- tibble(locus_rsid = region$rsid, trait1 = t1, trait2 = t2, n_snps = nrow(m))
  if (nrow(m) < cfg$coloc$min_snps) return(res)

  ds <- function(i, t) {
    n <- cfg$traits[[t]]$ncase + cfg$traits[[t]]$ncontrol
    list(snp = m$rsid, position = m$pos, beta = m[[paste0("beta", i)]],
         varbeta = m[[paste0("se", i)]]^2, type = "cc",
         N = n, s = cfg$traits[[t]]$ncase / n)
  }
  fit <- tryCatch(
    coloc.abf(ds(1, t1), ds(2, t2), p1 = cfg$coloc$p1, p2 = cfg$coloc$p2, p12 = cfg$coloc$p12),
    error = function(e) NULL)
  if (is.null(fit)) return(res)
  bind_cols(res, as_tibble(as.list(fit$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                                                  "PP.H3.abf", "PP.H4.abf")])))
}

res <- pmap_dfr(loci %>% select(rsid, traits), function(rsid, traits) {
  region <- regions[regions$rsid == rsid, ]
  if (length(traits) < 2) return(NULL)
  tests <- as_tibble(t(combn(sort(traits), 2)), .name_repair = ~ c("t1", "t2"))
  map2_dfr(tests$t1, tests$t2, coloc_pair, region = region)
}) %>%
  left_join(loci %>% select(locus_rsid = rsid, nearest_gene), by = "locus_rsid") %>%
  mutate(
    class = case_when(
      is.na(PP.H4.abf) ~ "Not tested",
      PP.H4.abf >= 0.8 ~ "Strong",
      PP.H4.abf >= 0.5 ~ "Moderate",
      PP.H3.abf > 0.5 ~ "Distinct",
      TRUE ~ "None"
    )
  )

write_tsv(res, file.path(out, "asset_coloc_abf.tsv"))

write_counts("coloc_abf", c(
  pairs_tested = sum(!is.na(res$PP.H4.abf)),
  pairs_h4_0.8 = sum(res$PP.H4.abf >= 0.8, na.rm = TRUE),
  loci_h4_0.8 = n_distinct(res$locus_rsid[res$PP.H4.abf >= 0.8 & !is.na(res$PP.H4.abf)])
))
