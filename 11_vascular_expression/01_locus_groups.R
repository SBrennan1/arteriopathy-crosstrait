# Loci and candidate genes for the cell-type analyses: ASSET pleiotropic loci, MTAG
# leads and trait-specific univariate leads (P >= 0.05 in the other traits), with every
# protein-coding gene within +/- 500 kb and the target gene flagged.

source("R/setup.R")
out <- res_dir("vascular")
p <- cfg$vascular
w <- p$locus_window_kb * 1e3
genes <- load_genes()

split_traits <- function(x) if (is.na(x)) character() else str_split_1(x, ",")

# Shared loci (ASSET) ----
asset <- read_stage("asset", "asset_loci.tsv")
v2g_asset <- read_stage("v2g", "asset_v2g.tsv") %>% select(locus_rsid, gene, evidence)

shared <- asset %>%
  left_join(v2g_asset, by = c(rsid = "locus_rsid")) %>%
  mutate(traits = map2_chr(traits_pos, traits_neg,
                           ~ paste(unique(c(split_traits(.x), split_traits(.y))), collapse = ","))) %>%
  transmute(locus = rsid, source = "ASSET", group = "Shared", chr, pos, traits, effect_pattern,
            target_gene = coalesce(gene, nearest_gene), gene_basis = coalesce(evidence, "Nearest gene"),
            nearest_gene)

# Trait-specific loci (univariate) ----
uni <- map(set_names(traits), read_univar, cols = c("rsid", "chr", "pos", "pval"))

specific <- map_dfr(traits, function(t) {
  leads <- clump_local(uni[[t]] %>% select(rsid, pval), p = cfg$traits[[t]]$gws_p)
  d <- uni[[t]] %>% filter(rsid %in% leads$rsid) %>% distinct(rsid, .keep_all = TRUE)
  others <- setdiff(traits, t)
  p_other <- map_dfc(set_names(others), function(o) {
    uni[[o]]$pval[match(d$rsid, uni[[o]]$rsid)]
  })
  d %>%
    mutate(n_other = rowSums(!is.na(p_other)),
           min_p_other = suppressWarnings(apply(p_other, 1, min, na.rm = TRUE)),
           trait = t)
}) %>%
  # Tested in at least three other traits and nominal in none
  filter(n_other >= 3, min_p_other >= p$specific_p_other)

near_shared <- function(chr, pos) {
  map2_lgl(chr, pos, ~ any(shared$chr == .x & abs(shared$pos - .y) <= w))
}
n_before <- nrow(specific)
specific <- specific %>%
  filter(!near_shared(chr, pos)) %>%
  bind_cols(nearest_gene(.$chr, .$pos, genes)) %>%
  transmute(locus = rsid, source = "Univariate", group = "Trait-specific", chr, pos,
            traits = trait, effect_pattern = "Single trait",
            target_gene = nearest_gene, gene_basis = "Nearest gene", nearest_gene)

# MTAG leads ----
mtag <- read_stage("mtag", "mtag_leads.tsv") %>%
  summarise(chr = first(chr), pos = first(pos), traits = paste(sort(unique(target_trait)), collapse = ","),
            .by = rsid) %>%
  inner_join(read_stage("v2g", "mtag_v2g.tsv") %>% select(rsid, gene, evidence, nearest_gene), by = "rsid") %>%
  transmute(locus = rsid, source = "MTAG", group = "MTAG", chr, pos, traits,
            effect_pattern = NA_character_,
            target_gene = coalesce(gene, nearest_gene), gene_basis = coalesce(evidence, "Nearest gene"),
            nearest_gene)

loci <- bind_rows(shared, specific, mtag)

# Genes in each locus window ----
locus_genes <- loci %>%
  select(locus, source, chr, pos, target_gene) %>%
  pmap_dfr(function(locus, source, chr, pos, target_gene) {
    genes %>%
      filter(chr == as.character(!!chr), end >= !!pos - w, start <= !!pos + w) %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      transmute(locus = !!locus, source = !!source, gene_id, gene = gene_name,
                is_target = gene_name == target_gene)
  })

missing_target <- loci %>% anti_join(filter(locus_genes, is_target), by = c("locus", "source"))
if (nrow(missing_target) > 0) {
  message(nrow(missing_target), " loci whose target gene is not a protein-coding gene in the window: ",
          paste(head(missing_target$target_gene, 10), collapse = ", "))
}

write_tsv(loci, file.path(out, "locus_groups.tsv"))
write_tsv(locus_genes, file.path(out, "locus_genes.tsv"))
write_tsv(distinct(locus_genes, gene_id, gene), file.path(out, "genes_for_expression.tsv"))

print(count(specific, traits))

write_counts("vascular_loci", c(
  shared = nrow(shared),
  specific_candidates = n_before,
  set_names(as.integer(table(factor(specific$traits, levels = traits))), paste0("specific_", traits)),
  mtag = nrow(mtag),
  genes = n_distinct(locus_genes$gene_id),
  targets_missing = nrow(missing_target)
))
