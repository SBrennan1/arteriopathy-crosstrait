# Target genes in aneurysm tissue: cell-type specificity of the IA target genes in the
# brain aneurysm atlas (Wang et al. 2026), separately in aneurysm and control donors, and
# differential expression of the AAD target genes in abdominal aortic aneurysm (GSE57691).
# Run 02_expression_summary.py brain_aneurysm first.

source("R/setup.R")
out <- res_dir("vascular")
p <- cfg$vascular
set.seed(p$seed)

loci <- read_stage("vascular", "locus_groups.tsv")
lg <- read_stage("vascular", "locus_genes.tsv")
score_units <- function(u, s) u %>% left_join(s, by = "gene") %>% select(unit, is_target, score)

# Brain aneurysm atlas ----
bd <- p$datasets$brain_aneurysm
cond <- p$brain_conditions
brain <- read_expr("brain_aneurysm") %>%
  mutate(assay = if_else(str_detect(donor, bd$sn_samples), "single-nucleus", "single-cell"))
donors <- distinct(brain, donor, region, assay)
print(count(donors, region, assay, name = "donors"))
write_tsv(donors, file.path(out, "brain_donors.tsv"))
u_ia <- units_for(lg, loci, c("ASSET", "MTAG", "Univariate"), trait = "IA")

# Per donor group: profile, and enrichment for cell types seen in at least 3 donors
brain_prof <- map_dfr(set_names(c(cond$case, cond$control)), function(cn) {
  celltype_profile(brain %>% filter(region == cn))
}, .id = "condition")
write_tsv(brain_prof, file.path(out, "brain_celltype_profile.tsv"))

brain_test <- function(units) {
  brain_prof %>%
    filter(n_donors >= 3) %>%
    distinct(condition, cell_type) %>%
    pmap_dfr(function(condition, cell_type) {
      # genes detectably expressed in this donor group (>= expressed_pct of one cell type)
      expressed <- brain_prof %>% filter(.data$condition == !!condition, pct >= p$expressed_pct) %>% distinct(gene)
      s <- brain_prof %>% filter(.data$condition == !!condition, .data$cell_type == !!cell_type) %>%
        semi_join(expressed, by = "gene") %>%
        select(gene, score = specificity)
      locus_perm(score_units(units, s)) %>% mutate(condition = condition, cell_type = cell_type, .before = 1)
    }) %>%
    add_q(condition)
}
brain_spec <- brain_test(u_ia)
write_tsv(brain_spec, file.path(out, "brain_ia_specificity.tsv"))
print(brain_spec %>% filter(q_greater < 0.05 | q_less < 0.05) %>%
        select(condition, cell_type, n_units, observed, null_mean, q_greater, q_less), n = Inf)

# Sensitivity: nearest-gene assignment
u_ia_nearest <- units_for(nearest_gene_targets(lg, loci), loci, c("ASSET", "MTAG", "Univariate"), trait = "IA")
brain_spec_nearest <- brain_test(u_ia_nearest)
write_tsv(brain_spec_nearest, file.path(out, "brain_ia_specificity_nearest_gene.tsv"))
print(brain_spec_nearest %>% filter(q_greater < 0.1 | q_less < 0.1) %>%
        select(condition, cell_type, n_units, observed, null_mean, q_greater, q_less), n = Inf)

# AAA bulk (GSE57691) ----
a <- p$aaa
dir.create(a$dir, recursive = TRUE, showWarnings = FALSE)
gse <- GEOquery::getGEO(a$geo, destdir = a$dir, GSEMatrix = TRUE)[[1]]
pd <- Biobase::pData(gse)
txt <- apply(pd[, grep("characteristics|source_name|title", names(pd)), drop = FALSE], 1, paste, collapse = " ")
grp <- case_when(str_detect(txt, regex(a$case_pattern, ignore_case = TRUE)) ~ "case",
                 str_detect(txt, regex(a$control_pattern, ignore_case = TRUE)) ~ "control")
print(table(grp, useNA = "ifany"))

x <- Biobase::exprs(gse)[, !is.na(grp)]
grp <- factor(grp[!is.na(grp)], levels = c("control", "case"))
if (max(x, na.rm = TRUE) > 50) x <- log2(pmax(x, 0) + 1)

fd <- Biobase::fData(gse)
sym_col <- intersect(c("Symbol", "ILMN_Gene", "Gene symbol", "GENE_SYMBOL"), names(fd))[1]
fit <- limma::eBayes(limma::lmFit(x, model.matrix(~ grp)))
aaa_de <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none") %>%
  as_tibble(rownames = "probe") %>%
  mutate(gene = fd[[sym_col]][match(probe, rownames(fd))]) %>%
  filter(!is.na(gene), gene != "") %>%
  slice_max(AveExpr, n = 1, with_ties = FALSE, by = gene) %>%
  transmute(gene, log_fc = logFC, t, p = P.Value, fdr = p.adjust(P.Value, "BH"), ave_expr = AveExpr)
write_tsv(aaa_de, file.path(out, "aaa_de.tsv"))

# AAD target genes
aa_targets <- loci %>%
  filter(map_lgl(str_split(traits, ","), ~ "AA" %in% .x)) %>%
  distinct(gene = target_gene, locus, source, group)
aaa_targets <- aa_targets %>%
  summarise(sources = paste(sort(unique(source)), collapse = ";"), .by = gene) %>%
  left_join(aaa_de, by = "gene") %>%
  arrange(fdr)
write_tsv(aaa_targets, file.path(out, "aaa_target_genes.tsv"))
print(aaa_targets %>% filter(fdr < 0.05) %>% select(gene, log_fc, fdr), n = Inf)

write_counts("vascular_aneurysm", c(
  brain_genes = n_distinct(brain$gene),
  brain_donors_aneurysm = sum(donors$region == cond$case),
  brain_donors_control = sum(donors$region == cond$control),
  ia_targets = n_distinct(u_ia$unit[u_ia$is_target]),
  aaa_cases = sum(grp == "case"),
  aaa_controls = sum(grp == "control"),
  aa_targets_on_array = sum(!is.na(aaa_targets$fdr)),
  aa_targets_fdr05 = sum(aaa_targets$fdr < 0.05, na.rm = TRUE)
))
