# Variant-to-gene mapping of the ASSET loci using the SuSiE credible sets (lead SNP
# where no credible set is available).

source("R/setup.R")
out <- res_dir("v2g")
p <- cfg$v2g

loci <- read_tsv(file.path(cfg$results_dir, "asset/asset_loci.tsv"), show_col_types = FALSE) %>%
  mutate(traits = map_chr(contributing_traits(traits_pos, traits_neg), paste, collapse = ","))
cs <- read_tsv(file.path(cfg$results_dir, "coloc/asset_credible_sets.tsv"), show_col_types = FALSE)
genes <- load_genes()

# A SNP in several credible sets keeps the PIP of the first trait listed
snps <- bind_rows(
  cs %>%
    distinct(locus_rsid, rsid, .keep_all = TRUE) %>%
    transmute(locus_rsid, rsid, chr, pos, chr_hg38, pos_hg38, a1 = effect_allele,
              a2 = other_allele, pip, source = "credible_set"),
  loci %>%
    filter(!rsid %in% cs$locus_rsid) %>%
    transmute(locus_rsid = rsid, rsid, chr, pos, chr_hg38, pos_hg38, a1 = A1, a2 = A2,
              pip = NA_real_, source = "lead_snp")
)

top <- snps %>%
  group_by(locus_rsid) %>%
  slice_max(coalesce(pip, 0), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(loci %>% select(locus_rsid = rsid, traits), by = "locus_rsid")

# Evidence tiers ----
vep <- vep_annotate(snps$rsid[coalesce(snps$pip, 0) > p$coding_min_pip])
coding <- snps %>%
  inner_join(vep, by = "rsid") %>%
  filter(is_coding, pip > p$coding_min_pip) %>%
  group_by(locus_rsid) %>%
  slice_max(pip, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(locus_rsid, coding_snp = rsid, coding_gene = vep_gene, coding_pip = pip,
            coding_consequence = vep_consequence)

coloc_all <- gtex_coloc(top %>% select(locus_rsid, chr, pos, traits))

e2g_snps <- snps %>% filter(pip > p$e2g_min_pip | source == "lead_snp")
e2g <- e2g_snps %>%
  inner_join(e2g_annotate(e2g_snps), by = "rsid") %>%
  filter(!is.na(e2g_gene), e2g_score > p$e2g_min_score) %>%
  group_by(locus_rsid) %>%
  slice_max(e2g_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(locus_rsid, e2g_snp = rsid, e2g_gene, e2g_score, e2g_tissue)

nearest <- top %>%
  transmute(locus_rsid, top_snp = rsid, top_pip = pip) %>%
  bind_cols(nearest_gene(top$chr, top$pos, genes))

# Assign genes ----
# gene_tiered: evidence tiers only; gene: after curated overrides (config/v2g_overrides.tsv)
v2g <- loci %>%
  transmute(locus_rsid = rsid, chr, pos, chr_hg38, pos_hg38, traits) %>%
  left_join(nearest, by = "locus_rsid") %>%
  left_join(count(filter(snps, source == "credible_set"), locus_rsid, name = "n_cs_snps"), by = "locus_rsid") %>%
  mutate(n_cs_snps = coalesce(n_cs_snps, 0L)) %>%
  left_join(coding, by = "locus_rsid") %>%
  left_join(best_eqtl_coloc(coloc_all, genes), by = "locus_rsid") %>%
  left_join(e2g, by = "locus_rsid") %>%
  assign_tiers() %>%
  apply_v2g_overrides("ASSET", key = "locus_rsid") %>%
  arrange(chr, pos)

write_tsv(v2g, file.path(out, "asset_v2g.tsv"))
write_tsv(coloc_all %>% mutate(gene = ensg_symbol(gene_id, genes), .after = gene_id),
          file.path(out, "asset_gtex_coloc.tsv"))

ev <- count(v2g, evidence)
write_counts("v2g_asset", c(
  loci = nrow(v2g),
  setNames(ev$n, paste0("evidence_", str_replace_all(tolower(ev$evidence), "[^a-z0-9]+", "_"))),
  functional_support = sum(v2g$evidence != "Nearest gene"),
  overridden = sum(v2g$overridden),
  eqtl_tests = nrow(coloc_all)
))
