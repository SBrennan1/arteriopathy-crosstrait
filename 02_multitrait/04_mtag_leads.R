# Pairwise MTAG lead variants (FUMA SNP2GENE) for each target trait, with univariate
# support and novelty relative to the univariate GWAS and to published associations
# (GWAS Catalog and config/known_loci_literature.tsv).

source("R/setup.R")
out <- res_dir("mtag")
fuma_in <- res_dir("mtag", "fuma_input")

pairs <- tibble(target = names(cfg$mtag$outputs), file = unlist(cfg$mtag$outputs)) %>%
  mutate(pair = sub("_trait_[12]\\.txt$", "", file),
         paired_trait = map2_chr(strsplit(pair, "_"), target, ~ setdiff(.x, .y)),
         pair = sub("_", "-", pair))

leads <- map_dfr(pairs$target, function(t) {
  fread(file.path(cfg$mtag$fuma_dir, t, "leadSNPs.txt"), colClasses = list(character = "rsID", numeric = "p")) %>%
    as_tibble() %>%
    transmute(rsid = rsID, chr = as.integer(chr), pos = as.integer(pos), p_fuma = p,
              target_trait = t)
}) %>%
  # Trait-specific significance threshold (stricter for AA)
  filter(p_fuma < map_dbl(target_trait, ~ cfg$traits[[.x]]$gws_p)) %>%
  select(-p_fuma) %>%
  left_join(pairs %>% select(target_trait = target, mtag_pair = pair, paired_trait),
            by = "target_trait")

# MTAG and univariate statistics at the lead SNPs, aligned to the univariate effect allele
lookup <- function(t, file = NULL) {
  d <- if (is.null(file)) read_univar(t) else vroom(file, show_col_types = FALSE)
  d %>% filter(rsid %in% leads$rsid) %>%
    select(rsid, effect_allele, other_allele, beta, pval) %>%
    mutate(trait = t)
}
univ <- map_dfr(traits, lookup)
mtag <- map_dfr(pairs$target, ~ lookup(.x, file.path(fuma_in, paste0(.x, "_mtag.txt.gz"))))

leads <- leads %>%
  left_join(univ %>% rename(target_trait = trait, effect_allele_univ = effect_allele,
                            other_allele_univ = other_allele, beta_univ = beta, p_univ = pval),
            by = c("rsid", "target_trait")) %>%
  left_join(mtag %>% rename(target_trait = trait, ea_mtag = effect_allele,
                            oa_mtag = other_allele, beta_mtag = beta, p_mtag = pval),
            by = c("rsid", "target_trait")) %>%
  left_join(univ %>% select(rsid, paired_trait = trait, p_paired_univ = pval),
            by = c("rsid", "paired_trait")) %>%
  mutate(
    flip = case_when(
      ea_mtag == effect_allele_univ & oa_mtag == other_allele_univ ~ 1,
      ea_mtag == other_allele_univ & oa_mtag == effect_allele_univ ~ -1
    ),
    beta_mtag = beta_mtag * flip,
    concordant = sign(beta_mtag) == sign(beta_univ),
    uplift = log10(p_univ) - log10(p_mtag),
    gws_univ = p_univ < map_dbl(target_trait, ~ cfg$traits[[.x]]$gws_p),
    novel_mtag = !gws_univ
  ) %>%
  select(-ea_mtag, -oa_mtag, -flip)

if (any(is.na(leads$concordant))) {
  warning(sum(is.na(leads$concordant)), " lead SNPs could not be allele-matched to the univariate GWAS")
}

# Univariate support: concordant direction (P threshold applied below)
supported <- leads %>% filter(p_univ < 0.05, concordant) %>% add_hg38()

# Previously reported associations for the same trait (GRCh38) ----
catalog <- fread(cfg$gwas_catalog, quote = "", select = c("DISEASE/TRAIT", "MAPPED_TRAIT",
                                                          "CHR_ID", "CHR_POS", "SNPS")) %>%
  as_tibble() %>%
  transmute(desc = paste(`DISEASE/TRAIT`, MAPPED_TRAIT), chr = CHR_ID,
            pos = suppressWarnings(as.integer(CHR_POS)), rsid = str_extract(SNPS, "rs[0-9]+")) %>%
  filter(!is.na(pos))

known <- map_dfr(traits, function(t) {
  catalog %>%
    filter(str_detect(desc, regex(paste(cfg$catalog_terms[[t]], collapse = "|"), ignore_case = TRUE))) %>%
    transmute(trait = t, chr, pos_hg38 = pos, rsid, source = "GWAS Catalog")
})

lit <- read_tsv(cfg$known_loci_extra, show_col_types = FALSE, col_types = cols(.default = "c"))
if (nrow(lit) > 0) {
  known <- bind_rows(known, lit %>%
                       mutate(pos = as.integer(pos)) %>%
                       add_hg38() %>%
                       transmute(trait, chr = chr_hg38, pos_hg38, rsid, source = reference))
}

window <- 500e3
novelty <- pmap_dfr(supported %>% select(rsid, target_trait, chr_hg38, pos_hg38),
                    function(rsid, target_trait, chr_hg38, pos_hg38) {
  k <- known %>% filter(trait == target_trait, chr == chr_hg38,
                        abs(pos_hg38 - !!pos_hg38) < 2 * window)
  near <- k %>% filter(abs(pos_hg38 - !!pos_hg38) < window)
  r2 <- ld_r2(c(rsid, k$rsid[!is.na(k$rsid)]))
  max_r2 <- if (!is.null(r2) && rsid %in% rownames(r2) && nrow(r2) > 1) {
    max(r2[rsid, setdiff(colnames(r2), rsid)])
  } else NA_real_
  tibble(rsid = rsid, target_trait = target_trait,
         known_within_500kb = nrow(near),
         known_min_dist_kb = if (nrow(k)) round(min(abs(k$pos_hg38 - pos_hg38)) / 1e3) else NA,
         known_max_r2 = max_r2,
         known_source = paste(unique(near$source), collapse = "; "))
})

supported <- supported %>%
  left_join(novelty, by = c("rsid", "target_trait")) %>%
  mutate(
    previously_reported = known_within_500kb > 0 | coalesce(known_max_r2 >= 0.1, FALSE),
    novel = novel_mtag & !previously_reported
  ) %>%
  arrange(target_trait, p_mtag)

write_tsv(leads, file.path(out, "mtag_leads_all.tsv"))
supported <- filter(supported, p_univ < cfg$mtag$univ_support_p)
write_tsv(supported, file.path(out, "mtag_leads.tsv"))

write_counts("mtag", c(
  lead_snp_trait_pairs = nrow(leads),
  supported_pairs = nrow(supported),
  supported_unique_snps = n_distinct(supported$rsid),
  not_gws_univariate_pairs = sum(supported$novel_mtag),
  novel_pairs = sum(supported$novel),
  novel_unique_snps = n_distinct(supported$rsid[supported$novel])
))
