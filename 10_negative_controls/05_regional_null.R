# Background rate of colocalization at the ASSET loci: each tested pair is retested with
# one member, chosen at random, replaced by its matched control in each of 100 draws
# (coloc.abf). Hair color is tested against every arteriopathy as a reference.

source("R/setup.R")
library(coloc)
nc <- cfg$negative_controls
out <- res_dir("negative_controls")
strong <- 0.8
set.seed(nc$seed)

draws <- fread(nc$draws)
n_draws <- nc$n_draws

# Control for trait t in draw d
control_for <- function(t, d) {
  ids <- draws[trait == t][order(draw)]$id
  ids[(d - 1) %% length(ids) + 1]
}
ctrl_n <- unique(draws[, .(id, ncase = n_cases_EUR, ncontrol = n_controls_EUR)])
setkey(ctrl_n, id)

read_control <- function(id) {
  f <- file.path(nc$dir, "regional", paste0(id, ".tsv"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, colClasses = list(character = c("chr", "ref", "alt")))
  d <- d[chr %in% as.character(1:22)][, chr := as.integer(chr)]
  setkey(d, chr, pos)
  d
}

# Arteriopathy data in the pleiotropic windows ----
w <- cfg$coloc$window_kb * 1e3
pleio <- read_tsv(file.path(cfg$results_dir, "asset/asset_loci.tsv"), show_col_types = FALSE)
win <- as.data.table(pleio)[, .(locus_rsid = rsid, chr = as.integer(chr),
                                start = as.integer(pmax(1, pos - w)), end = as.integer(pos + w))]

arterio <- map(set_names(traits), function(t) {
  d <- as.data.table(read_univar(t, cols = c("rsid", "chr", "pos", "effect_allele",
                                             "other_allele", "beta", "se")))
  d <- d[!is.na(beta) & !is.na(se) & se > 0 & !is.na(pos)]
  d <- d[win, on = .(chr, pos >= start, pos <= end), nomatch = 0,
         .(locus_rsid, rsid = x.rsid, chr = x.chr, pos = x.pos, effect_allele, other_allele, beta, se)]
  split(d, by = "locus_rsid", keep.by = FALSE)
})

cc <- function(beta, se, snp, ncase, ncontrol) {
  list(snp = snp, beta = beta, varbeta = se^2, type = "cc",
       N = ncase + ncontrol, s = ncase / (ncase + ncontrol))
}

coloc_null <- function(locus, kept, ctrl, id) {
  a <- arterio[[kept]][[locus]]
  if (is.null(a)) return(c(n_snps = 0, PP.H3.abf = NA, PP.H4.abf = NA))
  m <- ctrl[a, on = .(chr, pos), nomatch = 0]
  m <- m[(effect_allele == alt & other_allele == ref) | (effect_allele == ref & other_allele == alt)]
  m <- m[!duplicated(rsid)]
  if (nrow(m) < cfg$coloc$min_snps) return(c(n_snps = nrow(m), PP.H3.abf = NA, PP.H4.abf = NA))
  invisible(capture.output(fit <- suppressWarnings(suppressMessages(coloc.abf(
    cc(m$i.beta, m$i.se, m$rsid, cfg$traits[[kept]]$ncase, cfg$traits[[kept]]$ncontrol),
    cc(m$beta, m$se, m$rsid, ctrl_n[id]$ncase, ctrl_n[id]$ncontrol),
    p1 = cfg$coloc$p1, p2 = cfg$coloc$p2, p12 = cfg$coloc$p12)))))
  c(n_snps = nrow(m), fit$summary[c("PP.H3.abf", "PP.H4.abf")])
}

# Null coloc tests ----
obs <- read_tsv(file.path(cfg$results_dir, "coloc/asset_coloc_abf.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(PP.H4.abf)) %>%
  select(locus_rsid, trait1, trait2, PP.H4.abf)

tests <- expand_grid(obs %>% select(locus_rsid, trait1, trait2), draw = seq_len(n_draws)) %>%
  mutate(swap_first = runif(n()) < 0.5,
         kept = if_else(swap_first, trait2, trait1),
         swapped = if_else(swap_first, trait1, trait2),
         id = map2_chr(swapped, draw, control_for)) %>%
  select(-swap_first)

res <- map_dfr(unique(tests$id), function(id) {
  ctrl <- read_control(id)
  if (is.null(ctrl)) {
    warning("no regional file for ", id, ", rerun 02_extract_regions.sh", call. = FALSE)
    return(NULL)
  }
  g <- tests[tests$id == id, ]
  bind_cols(g, map_dfr(seq_len(nrow(g)), ~ as_tibble(as.list(coloc_null(g$locus_rsid[.x], g$kept[.x], ctrl, id)))))
})
fwrite(res, file.path(out, "coloc_null_tests.tsv.gz"), sep = "\t")

coloc_by_draw <- res %>%
  group_by(draw) %>%
  summarise(pairs_tested = sum(!is.na(PP.H4.abf)),
            strong_pairs = sum(PP.H4.abf >= strong, na.rm = TRUE),
            strong_loci = n_distinct(locus_rsid[PP.H4.abf >= strong & !is.na(PP.H4.abf)]))
write_tsv(coloc_by_draw, file.path(out, "coloc_null_by_draw.tsv"))

null_summary <- function(observed, null) {
  tibble(observed = observed, null_mean = mean(null),
         null_lower = quantile(null, 0.025), null_upper = quantile(null, 0.975),
         empirical_p = (1 + sum(null >= observed)) / (length(null) + 1))
}
coloc_summary <- bind_rows(
  null_summary(sum(obs$PP.H4.abf >= strong), coloc_by_draw$strong_pairs) %>% mutate(measure = "strong_pairs"),
  null_summary(n_distinct(obs$locus_rsid[obs$PP.H4.abf >= strong]), coloc_by_draw$strong_loci) %>%
    mutate(measure = "strong_loci"),
  null_summary(mean(obs$PP.H4.abf >= strong), coloc_by_draw$strong_pairs / coloc_by_draw$pairs_tested) %>%
    mutate(measure = "strong_pair_rate")
) %>%
  mutate(pairs_observed = nrow(obs), .before = 1) %>%
  relocate(measure)
write_tsv(coloc_summary, file.path(out, "coloc_null_summary.tsv"))
print(coloc_summary)

# Rate by arteriopathy
coloc_by_trait <- map_dfr(traits, function(t) {
  o <- obs %>% filter(trait1 == t | trait2 == t)
  n <- res %>% filter(swapped == t, !is.na(PP.H4.abf))
  tibble(trait = t, pairs = nrow(o), observed_rate = mean(o$PP.H4.abf >= strong),
         null_tests = nrow(n), null_rate = mean(n$PP.H4.abf >= strong))
})
write_tsv(coloc_by_trait, file.path(out, "coloc_null_by_trait.tsv"))

# Reference trait (hair color) ----
ref <- if (file.exists(nc$reference_file)) fread(nc$reference_file) else data.table()
ref_sum <- tibble()
if (nrow(ref) > 0) {
  ctrl_n <- unique(rbind(ctrl_n, ref[, .(id, ncase = n_cases_EUR, ncontrol = n_controls_EUR)]))
  setkey(ctrl_n, id)
  ref_tests <- bind_rows(obs %>% select(locus_rsid, trait = trait1),
                         obs %>% select(locus_rsid, trait = trait2)) %>%
    distinct()
  ref_res <- map_dfr(seq_len(nrow(ref)), function(i) {
    ctrl <- read_control(ref$id[i])
    if (is.null(ctrl)) {
      warning("no regional file for ", ref$id[i], ", rerun 02_extract_regions.sh", call. = FALSE)
      return(NULL)
    }
    bind_cols(ref_tests %>% mutate(label = ref$label[i], .before = 1),
              map_dfr(seq_len(nrow(ref_tests)), ~ as_tibble(as.list(
                coloc_null(ref_tests$locus_rsid[.x], ref_tests$trait[.x], ctrl, ref$id[i])))))
  })
  if (nrow(ref_res) > 0) {
    write_tsv(ref_res, file.path(out, "coloc_reference.tsv"))
    ref_sum <- ref_res %>%
      group_by(label) %>%
      summarise(tests = sum(!is.na(PP.H4.abf)),
                n_strong = sum(PP.H4.abf >= strong, na.rm = TRUE),
                n_strong_loci = n_distinct(locus_rsid[PP.H4.abf >= strong & !is.na(PP.H4.abf)]),
                max_pph4 = suppressWarnings(max(PP.H4.abf, na.rm = TRUE)))
    print(ref_sum)
  }
}

write_counts("negative_controls_regional", c(
  coloc_strong_pairs_observed = coloc_summary$observed[1],
  coloc_strong_pairs_null_mean = round(coloc_summary$null_mean[1], 2),
  coloc_strong_loci_observed = coloc_summary$observed[2],
  coloc_strong_loci_null_mean = round(coloc_summary$null_mean[2], 2),
  coloc_empirical_p = signif(coloc_summary$empirical_p[1], 3),
  if (nrow(ref_sum)) set_names(c(ref_sum$tests, ref_sum$n_strong),
                               c(paste0("reference_tests_", ref_sum$label), paste0("reference_strong_", ref_sum$label)))
))
