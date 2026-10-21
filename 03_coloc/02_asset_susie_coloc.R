# SuSiE fine-mapping (UK Biobank LD, +/- 500 kb) of each trait contributing to an
# ASSET pleiotropic locus, and coloc.susie for each pair of contributing traits
# (sensitivity analysis).

source("R/setup.R")
out <- res_dir("coloc")
p <- cfg$susie

ukbb_ld_init()
catalogue <- ukbb_ld_catalogue()

w <- p$window_kb * 1e3
loci <- read_tsv(file.path(cfg$results_dir, "asset/asset_loci.tsv"), show_col_types = FALSE) %>%
  mutate(traits = contributing_traits(traits_pos, traits_neg),
         start = pmax(1, pos - w), end = pos + w)

# Load data ----
# Median N of each GWAS used for susie_rss
windows <- data.table(chr = as.integer(loci$chr), start = as.integer(loci$start),
                      end = as.integer(loci$end))
sumstats <- map(set_names(unique(unlist(loci$traits))), function(t) {
  d <- read_univar(t, cols = c("rsid", "chr", "pos", "effect_allele", "other_allele",
                               "beta", "se", "pval", "N"))
  n <- median(d$N, na.rm = TRUE)
  d <- distinct(d, rsid, .keep_all = TRUE)
  hit <- as.data.table(d)[windows, on = .(chr, pos >= start, pos <= end), nomatch = 0, which = TRUE]
  list(d = d[sort(unique(hit)), ], n = n)
})

# Fine-mapping and coloc.susie per locus ----
run_locus <- function(i) {
  l <- loci[i, ]
  tr <- l$traits[[1]]
  pairs <- if (length(tr) >= 2) combn(sort(tr), 2, simplify = FALSE) else list()
  need_abf <- function(t1, t2, reason) tibble(locus_rsid = l$rsid, trait1 = t1, trait2 = t2, reason = reason)

  ld <- tryCatch(load_ukbb_ld(catalogue, l$chr, l$start, l$end), error = function(e) NULL)
  if (is.null(ld)) {
    return(list(
      finemap = tibble(locus_rsid = l$rsid, trait = tr, status = "no UKBB LD chunk"),
      abf = map_dfr(pairs, ~ need_abf(.x[1], .x[2], "no UKBB LD chunk"))
    ))
  }

  fits <- map(set_names(tr), function(t) {
    region <- sumstats[[t]]$d %>% filter(chr == l$chr, pos >= l$start, pos <= l$end)
    finemap_trait(region, ld, sumstats[[t]]$n)
  })

  coloc <- list()
  abf <- list()
  for (pr in pairs) {
    f1 <- fits[[pr[1]]]
    f2 <- fits[[pr[2]]]
    why <- case_when(
      is.null(f1$fit) ~ paste0(pr[1], ": ", f1$status),
      is.null(f2$fit) ~ paste0(pr[2], ": ", f2$status),
      !f1$has_cs ~ paste0(pr[1], ": no credible set"),
      !f2$has_cs ~ paste0(pr[2], ": no credible set"),
      TRUE ~ NA_character_
    )
    co <- if (is.na(why)) coloc_susie_pair(f1, f2)
    if (is.null(co)) {
      abf[[length(abf) + 1]] <- need_abf(pr[1], pr[2], coalesce(why, "coloc.susie failed"))
      next
    }
    coloc[[length(coloc) + 1]] <- tibble(
      locus_rsid = l$rsid, trait1 = pr[1], trait2 = pr[2],
      n_snps = min(f1$n_in_ld, f2$n_in_ld), co,
      lead_snp_t1 = f1$lead_snp, lead_p_t1 = f1$lead_p,
      lead_snp_t2 = f2$lead_snp, lead_p_t2 = f2$lead_p,
      cs_size_t1 = f1$cs_size, cs_size_t2 = f2$cs_size
    )
  }

  list(
    finemap = map_dfr(fits, finemap_summary, .id = "trait") %>% mutate(locus_rsid = l$rsid, .before = 1),
    cs = map_dfr(fits, "credible", .id = "trait") %>% mutate(locus_rsid = l$rsid, .before = 1),
    coloc = bind_rows(coloc),
    abf = bind_rows(abf)
  )
}

res <- map(seq_len(nrow(loci)), run_locus)

locus_pos <- loci %>% select(locus_rsid = rsid, nearest_gene, chr, pos, chr_hg38, pos_hg38)

coloc_tbl <- map_dfr(res, "coloc") %>%
  mutate(class = coloc_class(PP.H4.abf, PP.H3.abf), method = "coloc.susie", ld_source = "UKBB") %>%
  left_join(locus_pos, by = "locus_rsid") %>%
  relocate(nearest_gene, chr, pos, chr_hg38, pos_hg38, .after = locus_rsid) %>%
  arrange(desc(PP.H4.abf))

cs_tbl <- map_dfr(res, "cs") %>% add_hg38()
finemap_tbl <- map_dfr(res, "finemap")

abf_tbl <- map_dfr(res, "abf") %>%
  left_join(loci %>% select(locus_rsid = rsid, chr, pos, chr_hg38, pos_hg38, start, end), by = "locus_rsid") %>%
  relocate(chr, pos, chr_hg38, pos_hg38, start, end, .after = locus_rsid)

write_tsv(coloc_tbl, file.path(out, "asset_susie_coloc.tsv"))
write_tsv(cs_tbl, file.path(out, "asset_credible_sets.tsv"))
write_tsv(finemap_tbl, file.path(out, "asset_susie_finemap.tsv"))
write_tsv(abf_tbl, file.path(out, "asset_coloc_need_abf.tsv"))

write_counts("coloc_susie", c(
  loci = nrow(loci),
  loci_with_ld = n_distinct(finemap_tbl$locus_rsid[finemap_tbl$status != "no UKBB LD chunk"]),
  trait_fits = sum(finemap_tbl$status == "ok"),
  trait_fits_with_cs = sum(finemap_tbl$n_cs > 0, na.rm = TRUE),
  pairs_coloc_susie = nrow(coloc_tbl),
  pairs_h4_0.8 = sum(coloc_tbl$PP.H4.abf >= p$pph4_strong, na.rm = TRUE),
  loci_h4_0.8 = n_distinct(coloc_tbl$locus_rsid[coloc_tbl$PP.H4.abf >= p$pph4_strong & !is.na(coloc_tbl$PP.H4.abf)]),
  pairs_need_abf = nrow(abf_tbl)
))
