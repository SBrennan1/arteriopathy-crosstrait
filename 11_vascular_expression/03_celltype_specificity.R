# Cell-type specificity of the target genes of the ASSET and MTAG loci in the HCA arterial
# atlas (Zhao et al. 2025), compared with other expressed genes at the same loci.
# Sensitivity analyses: expression-matched comparison genes; nearest-gene assignment.
# Run 02_expression_summary.py hca_artery first.

source("R/setup.R")
out <- res_dir("vascular")
p <- cfg$vascular
set.seed(p$seed)

loci <- read_stage("vascular", "locus_groups.tsv")
lg <- read_stage("vascular", "locus_genes.tsv")
expr <- read_expr("hca_artery")

prof <- celltype_profile(expr)
cell_types <- sort(unique(prof$cell_type))
write_tsv(prof, file.path(out, "hca_celltype_profile.tsv"))

# Target genes: one row per gene with its sources and specificity per cell type
targets <- lg %>%
  filter(is_target) %>%
  left_join(loci %>% select(locus, source, group, traits), by = c("locus", "source")) %>%
  summarise(sources = paste(sort(unique(source)), collapse = ";"),
            traits = paste(sort(unique(unlist(str_split(traits, ",")))), collapse = ";"),
            .by = c(gene_id, gene))

target_spec <- prof %>%
  filter(gene_id %in% targets$gene_id) %>%
  mutate(expressed = pct >= p$expressed_pct) %>%
  summarise(top_cell_type = cell_type[which.max(specificity)],
            top_specificity = max(specificity),
            tau = first(tau),
            max_pct = max(pct),
            expressed_in = paste(cell_type[expressed], collapse = ";"),
            .by = c(gene_id, gene)) %>%
  left_join(prof %>%
              filter(gene_id %in% targets$gene_id) %>%
              select(gene_id, cell_type, specificity) %>%
              pivot_wider(names_from = cell_type, values_from = specificity, names_prefix = "spec_"),
            by = "gene_id") %>%
  right_join(targets, by = c("gene_id", "gene")) %>%
  arrange(desc(top_specificity))
write_tsv(target_spec, file.path(out, "hca_target_specificity.tsv"))

# Enrichment against the other genes in the same loci ----
# Mean specificity of the target genes against 10,000 draws of one protein-coding gene
# from each target's locus window; genes expressed in >= 10% of cells of a cell type.
score_units <- function(u, score_tbl) {
  u %>% left_join(score_tbl, by = "gene") %>% select(unit, is_target, score)
}
u <- units_for(lg, loci, c("ASSET", "MTAG"))
expressed_any <- prof %>%
  summarise(score = as.numeric(any(pct >= p$expressed_pct)), .by = gene)
expressed_genes <- expressed_any$gene[expressed_any$score == 1]
spec_of <- function(ct) prof %>% filter(cell_type == ct, gene %in% expressed_genes) %>% select(gene, score = specificity)
by_cell_type <- function(units) {
  map_dfr(set_names(cell_types), ~ locus_perm(score_units(units, spec_of(.x))), .id = "cell_type") %>% add_q()
}

# Expressed in any cell type
expressed_test <- locus_perm(score_units(u, expressed_any)) %>% mutate(cell_type = "Expressed (any cell type)")
enrich <- bind_rows(by_cell_type(u), expressed_test)
write_tsv(enrich, file.path(out, "hca_celltype_enrichment.tsv"))

# Sensitivity: comparison genes matched on expression level
level <- prof %>% summarise(level = mean(mean_log), .by = gene)
genes_main <- u %>% summarise(is_target = any(is_target), .by = gene) %>% left_join(level, by = "gene")
enrich_expr <- map_dfr(set_names(cell_types), function(ct) {
  expr_matched_perm(genes_main %>% left_join(spec_of(ct), by = "gene"))
}, .id = "cell_type") %>%
  add_q()
write_tsv(enrich_expr, file.path(out, "hca_celltype_enrichment_exprmatched.tsv"))

# Sensitivity: nearest-gene assignment
eqtl_loci <- loci %>% filter(source %in% c("ASSET", "MTAG"), gene_basis == "eQTL colocalization")
message(nrow(eqtl_loci), " loci assigned by eQTL colocalization, of which ",
        sum(eqtl_loci$target_gene == eqtl_loci$nearest_gene, na.rm = TRUE), " to the nearest gene")
u_nearest <- units_for(nearest_gene_targets(lg, loci), loci, c("ASSET", "MTAG"))
enrich_nearest <- by_cell_type(u_nearest)
write_tsv(enrich_nearest, file.path(out, "hca_celltype_enrichment_nearest_gene.tsv"))

# Dot plot data (Figure 5a) ----
write_tsv(prof %>% filter(gene_id %in% targets$gene_id), file.path(out, "hca_dotplot_celltype.tsv"))

print(enrich %>% select(cell_type, n_units, observed, null_mean, p_greater, q_greater, p_less, q_less))
print(enrich_expr %>% select(cell_type, observed, null_mean, p_greater, q_greater, q_less))
print(enrich_nearest %>% select(cell_type, n_units, observed, null_mean, p_greater, q_greater, q_less))

write_counts("vascular_hca_celltype", c(
  genes_in_atlas = n_distinct(prof$gene_id),
  targets = nrow(targets),
  targets_in_atlas = sum(!is.na(target_spec$top_specificity)),
  targets_expressed = sum(coalesce(target_spec$expressed_in, "") != "")
))
