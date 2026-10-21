# Genotype-Phenotype Map (GPMap) colocalization of the five arteriopathy GWAS, MR of
# blood pressure on each arteriopathy, and direction of effect relative to CAD
# (Aragam et al. 2022). Figure 3. Requires OPENGWAS_JWT for the blood pressure GWAS.

source("R/setup.R")
library(ggplot2)
library(cowplot)
library(ggrepel)
library(TwoSampleMR)
out <- res_dir("gpmap")
p <- cfg$gpmap
art_order <- c("AA", "IA", "SCAD", "FMD", "CeAD")

if (!nzchar(Sys.getenv("OPENGWAS_JWT"))) stop("Set OPENGWAS_JWT to an OpenGWAS token")

# GPMap results ----
# Results are cached, as GPMap is updated over time
cg_file <- file.path(out, "gpmap_coloc_groups.rds")
if (!file.exists(cg_file)) {
  imap(p$guids, ~ gpmapr::trait(.x)$coloc_groups %>% mutate(arteriopathy = .y)) %>%
    bind_rows() %>%
    saveRDS(cg_file)
}
all_cg <- as_tibble(readRDS(cg_file))

simplify_trait <- function(name) {
  name <- gsub("\\(UKB data field \\d+\\)", "", name)
  name <- gsub("\\(Firth correction\\)|\\(SPA correction\\)", "", name)
  name <- gsub(", automated reading|automated reading", "", name)
  name <- gsub("Diagnoses - (main|secondary) ICD10: [A-Z]\\d+\\.?\\d* ", "", name)
  name <- gsub("Non-cancer illness code, self-reported: ", "", name)
  name <- gsub("Vascular/heart problems diagnosed by doctor: ", "", name)
  name <- gsub("Operative procedures.*OPCS: [A-Z]\\d+\\.\\d+ ", "", name)
  name <- gsub("Medication for pain relief, constipation, heartburn: ", "Medication: ", name)
  name <- gsub("Treatment/medication code: ", "Medication: ", name)
  name <- gsub("Illnesses of (father|mother|siblings): ", "Family history: ", name)
  name <- gsub("Medication use \\(", "Medication: ", name)
  name <- gsub("\\)$", "", name)
  gsub("\\s+", " ", trimws(name))
}

# Complex traits (Phenotype) and molecular QTLs
complex <- all_cg %>%
  filter(data_type == "Phenotype") %>%
  mutate(trait_simple = simplify_trait(trait_name),
         trait_category = if_else(trait_simple == "Atherosclerotic heart disease" &
                                    trait_category == "Disease Of Ear And Mastoid Process",
                                  "Disease Of Circulatory System", trait_category))
molecular <- all_cg %>% filter(!is.na(data_type), data_type != "Phenotype")

# Locus annotation ----
# Gene from GPMap molecular colocalization (expression > protein > splicing >
# methylation), otherwise the nearest protein-coding gene
mol_priority <- c("Gene Expression" = 1, "Protein" = 2, "Splice Variant" = 3, "Methylation" = 4)
mol_annot <- molecular %>%
  filter(!is.na(gene), gene != "") %>%
  mutate(priority = mol_priority[data_type]) %>%
  group_by(arteriopathy, coloc_group_id) %>%
  summarise(best_gene = gene[order(priority)][1],
            best_tissue = tissue[order(priority)][1],
            best_data_type = data_type[order(priority)][1],
            all_genes = paste(unique(gene), collapse = "; "),
            all_tissues = paste(unique(tissue), collapse = "; "),
            all_data_types = paste(unique(data_type), collapse = "; "),
            n_molecular_colocs = n(), .groups = "drop")

complex_summary <- complex %>%
  group_by(arteriopathy, coloc_group_id) %>%
  summarise(n_complex_traits = n_distinct(trait_name),
            trait_categories = paste(unique(trait_category[!is.na(trait_category)]), collapse = "; "),
            top_traits = paste(head(unique(trait_name), 5), collapse = "; "), .groups = "drop")

# One row per arteriopathy colocalization group
locus_table <- all_cg %>%
  filter(is.na(data_type)) %>%
  distinct(arteriopathy, coloc_group_id, display_snp, ld_block_id, chr, bp) %>%
  left_join(mol_annot, by = c("arteriopathy", "coloc_group_id")) %>%
  left_join(complex_summary, by = c("arteriopathy", "coloc_group_id")) %>%
  mutate(has_molecular = !is.na(best_gene),
         snp_chr = as.numeric(sub(":.*", "", display_snp)),
         snp_pos = as.numeric(sub(".*:(\\d+)\\s.*", "\\1", display_snp)))

unannotated <- locus_table %>%
  filter(!has_molecular, !is.na(snp_pos)) %>%
  distinct(display_snp, snp_chr, snp_pos)
nearest <- unannotated %>%
  select(display_snp) %>%
  bind_cols(nearest_gene(unannotated$snp_chr, unannotated$snp_pos, load_genes()))

# Gene labels for lncRNA, pseudogene and Ensembl IDs (config: gpmap$label_fixes)
label_fixes <- unlist(p$label_fixes)
locus_table <- locus_table %>%
  left_join(nearest, by = "display_snp") %>%
  mutate(mapped_gene = if_else(has_molecular, best_gene, nearest_gene),
         mapping_source = if_else(has_molecular, paste0("GPMap molecular (", best_data_type, ")"),
                                  "Nearest protein-coding gene"),
         mapped_gene = coalesce(unname(label_fixes[mapped_gene]), mapped_gene))

# Variant-to-gene assignment of the matching ASSET locus
asset_by_label <- read_stage("asset", "asset_loci.tsv") %>%
  select(locus_rsid = rsid, asset_label = nearest_gene) %>%
  inner_join(read_stage("v2g", "asset_v2g.tsv") %>%
               select(locus_rsid, asset_v2g_gene = gene, asset_v2g_evidence = evidence),
             by = "locus_rsid") %>%
  distinct(asset_label, .keep_all = TRUE) %>%
  select(-locus_rsid)
locus_table <- left_join(locus_table, asset_by_label, by = c(mapped_gene = "asset_label"))

first_label <- locus_table %>% distinct(display_snp, .keep_all = TRUE)
gene_map <- setNames(first_label$mapped_gene, first_label$display_snp)
gene_label <- function(snp) if_else(snp %in% names(gene_map), unname(gene_map[snp]), snp)

# Supplementary tables ----
bp_pattern <- "blood pressure|hypertension|systolic|diastolic|pulse pressure"

locus_table %>%
  transmute(arteriopathy, candidate_variant = display_snp, chr = snp_chr, pos = snp_pos,
            mapped_gene, mapping_source, gpmap_all_genes = all_genes, gpmap_tissues = all_tissues,
            gpmap_data_types = all_data_types, n_molecular_colocs, n_complex_traits,
            trait_categories, top_colocalizing_traits = top_traits, nearest_gene,
            nearest_distance, asset_v2g_gene, asset_v2g_evidence) %>%
  arrange(arteriopathy, chr, pos) %>%
  write_tsv(file.path(out, "gpmap_locus_annotations.tsv"))

complex %>%
  filter(!is.na(trait_category)) %>%
  group_by(trait_category, arteriopathy) %>%
  summarise(n_loci = n_distinct(coloc_group_id), .groups = "drop") %>%
  pivot_wider(names_from = arteriopathy, values_from = n_loci, values_fill = 0) %>%
  select(trait_category, any_of(art_order)) %>%
  mutate(total = rowSums(across(any_of(art_order)))) %>%
  arrange(desc(total)) %>%
  write_tsv(file.path(out, "gpmap_trait_category_matrix.tsv"))

complex %>%
  group_by(trait_name, trait_category) %>%
  summarise(n_arteriopathies = n_distinct(arteriopathy),
            which_arteriopathies = paste(sort(unique(arteriopathy)), collapse = ", "),
            n_total_loci = n_distinct(paste(arteriopathy, coloc_group_id)), .groups = "drop") %>%
  filter(n_arteriopathies >= 3) %>%
  arrange(desc(n_arteriopathies), desc(n_total_loci)) %>%
  write_tsv(file.path(out, "gpmap_traits_3plus.tsv"))

# Blood pressure colocalizations grouped by locus gene
complex %>%
  filter(grepl(bp_pattern, trait_name, ignore.case = TRUE)) %>%
  mutate(gene = gene_label(display_snp)) %>%
  group_by(gene) %>%
  summarise(candidate_variants = paste(unique(display_snp), collapse = "; "),
            arteriopathies = paste(sort(unique(arteriopathy)), collapse = ", "),
            n_arteriopathies = n_distinct(arteriopathy),
            bp_traits = paste(unique(trait_name), collapse = "; "), .groups = "drop") %>%
  arrange(desc(n_arteriopathies)) %>%
  select(gene, candidate_variants, n_arteriopathies, arteriopathies, bp_traits) %>%
  write_tsv(file.path(out, "gpmap_bp_loci.tsv"))

molecular %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarise(n_arteriopathies = n_distinct(arteriopathy),
            which = paste(sort(unique(arteriopathy)), collapse = ", "),
            tissues = paste(unique(tissue), collapse = "; "),
            data_types = paste(unique(data_type), collapse = "; "), .groups = "drop") %>%
  arrange(desc(n_arteriopathies)) %>%
  write_tsv(file.path(out, "gpmap_gene_molecular_overlap.tsv"))

# MR: blood pressure -> arteriopathies ----
bp_ids <- unlist(p$bp_mr$ids)
bp_info <- as_tibble(ieugwasr::gwasinfo(unname(bp_ids)))
write_tsv(bp_info %>% select(any_of(c("id", "trait", "unit", "sample_size", "nsnp", "author",
                                      "consortium", "year", "pmid", "population", "build"))),
          file.path(out, "bp_exposures.tsv"))

units <- tolower(bp_info$unit[match(bp_ids[c("SBP", "DBP")], bp_info$id)])
if (all(str_detect(coalesce(units, ""), "mm ?hg"))) {
  bp_scale <- p$bp_mr$or_per_mmhg
  bp_unit <- "mmHg"
  x_lab <- sprintf("OR per %g mmHg increase in BP (95%% CI)", bp_scale)
} else if (all(str_detect(coalesce(units, ""), "\\bsd\\b|standard deviation"))) {
  bp_scale <- 1
  bp_unit <- "SD"
  x_lab <- "OR per SD increase in BP (95% CI)"
} else {
  stop("Unrecognised units for the BP exposure GWAS: ", paste(units, collapse = ", "))
}

instruments <- imap_dfr(bp_ids, function(id, bp) {
  extract_instruments(id, p1 = p$bp_mr$p, clump = TRUE, r2 = p$bp_mr$clump_r2,
                      kb = p$bp_mr$clump_kb) %>%
    mutate(bp_trait = bp)
})

outcomes <- map(set_names(art_order), function(t) {
  read_mr_outcome(t, unique(instruments$SNP)) %>%
    as.data.frame() %>%
    format_data(type = "outcome", snp_col = "rsid") %>%
    mutate(outcome = t, id.outcome = t)
})

bp_mr <- function(bp, t) {
  harm <- harmonise_data(filter(instruments, bp_trait == bp), outcomes[[t]]) %>% filter(mr_keep)
  if (nrow(harm) < p$bp_mr$min_snps) return(NULL)
  tag <- function(d) as_tibble(d) %>% mutate(bp_trait = bp, arteriopathy = t)
  list(mr = tag(mr(harm, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                                         "mr_weighted_mode"))),
       het = tag(mr_heterogeneity(harm)),
       pleio = tag(mr_pleiotropy_test(harm)))
}
mr_res <- expand_grid(bp = names(bp_ids), t = art_order) %>%
  pmap(bp_mr) %>%
  compact()

mr_all <- map_dfr(mr_res, "mr")
mr_all %>%
  mutate(exposure_unit = bp_unit, or = exp(b), or_lci = exp(b - 1.96 * se), or_uci = exp(b + 1.96 * se)) %>%
  select(bp_trait, arteriopathy, method, nsnp, b, se, exposure_unit, or, or_lci, or_uci, pval) %>%
  arrange(bp_trait, arteriopathy, method) %>%
  write_tsv(file.path(out, "bp_mr_all_methods.tsv"))
map_dfr(mr_res, "het") %>%
  select(bp_trait, arteriopathy, method, Q, Q_df, Q_pval) %>%
  write_tsv(file.path(out, "bp_mr_heterogeneity.tsv"))
map_dfr(mr_res, "pleio") %>%
  select(bp_trait, arteriopathy, egger_intercept, se, pval) %>%
  write_tsv(file.path(out, "bp_mr_egger_intercept.tsv"))

# Direction of effect relative to CAD ----
# GPMap lead variants aligned to the CAD risk-increasing allele; classified where both
# are genome-wide significant
coloc_variants <- complex %>%
  filter(!is.na(rsid), rsid != "") %>%
  distinct(arteriopathy, coloc_group_id, ld_block_id, display_snp, rsid, min_p) %>%
  left_join(locus_table %>% distinct(display_snp, mapped_gene), by = "display_snp") %>%
  mutate(locus_label = if_else(!is.na(mapped_gene) & mapped_gene != "", mapped_gene, display_snp))

art_effects <- map_dfr(art_order, function(t) {
  read_univar(t, cols = c("rsid", "chr", "pos", "effect_allele", "other_allele", "beta", "se", "pval")) %>%
    filter(rsid %in% coloc_variants$rsid) %>%
    group_by(rsid) %>%
    slice_min(pval, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(rsid, chr, pos, arteriopathy = t, art_beta = beta, art_se = se, art_pval = pval,
              art_ea = effect_allele, art_oa = other_allele)
})

# GRCh37 positions for the CAD lookup
pos37 <- art_effects %>%
  arrange(arteriopathy != "AA") %>%
  distinct(rsid, .keep_all = TRUE) %>%
  select(rsid, chr, pos)

cad <- vroom(p$cad_sumstats, show_col_types = FALSE,
             col_select = c(chromosome, base_pair_location, effect_allele, other_allele,
                            beta, standard_error, p_value)) %>%
  filter(as.character(chromosome) %in% as.character(1:22)) %>%
  mutate(chromosome = as.integer(chromosome), base_pair_location = as.integer(base_pair_location)) %>%
  inner_join(pos37, by = c(chromosome = "chr", base_pair_location = "pos")) %>%
  transmute(rsid, cad_beta = as.numeric(beta), cad_se = as.numeric(standard_error),
            cad_pval = as.numeric(p_value), cad_ea = toupper(effect_allele),
            cad_oa = toupper(other_allele)) %>%
  mutate(cad_risk_allele = if_else(cad_beta > 0, cad_ea, cad_oa),
         cad_nonrisk_allele = if_else(cad_beta > 0, cad_oa, cad_ea))

strand_flip <- function(a) chartr("ACGT", "TGCA", a)
palindromic <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
gws <- 5e-8

dir_pairs <- coloc_variants %>%
  left_join(art_effects %>% select(-chr, -pos), by = c("rsid", "arteriopathy")) %>%
  left_join(cad, by = "rsid") %>%
  filter(!is.na(art_beta), !is.na(cad_beta)) %>%
  mutate(
    art_beta_aligned = case_when(
      art_ea == cad_risk_allele ~ art_beta,
      art_ea == cad_nonrisk_allele ~ -art_beta,
      strand_flip(art_ea) == cad_risk_allele ~ art_beta,
      strand_flip(art_ea) == cad_nonrisk_allele ~ -art_beta
    ),
    art_z_aligned = art_beta_aligned / art_se,
    alignment_status = case_when(
      art_ea == cad_risk_allele | art_ea == cad_nonrisk_allele ~ "direct",
      strand_flip(art_ea) == cad_risk_allele | strand_flip(art_ea) == cad_nonrisk_allele ~ "strand_flip",
      TRUE ~ "unresolved"
    ),
    needs_manual_check = palindromic(art_ea, art_oa) & alignment_status == "strand_flip",
    direction = case_when(
      is.na(art_beta_aligned) ~ "unresolved",
      !(art_pval < gws & cad_pval < gws) ~ "not_both_sig",
      art_beta_aligned > 0 ~ "concordant",
      art_beta_aligned < 0 ~ "discordant",
      TRUE ~ "unresolved"
    )
  )

# One point per locus and arteriopathy
directionality <- dir_pairs %>%
  filter(direction %in% c("concordant", "discordant")) %>%
  group_by(locus_label, arteriopathy) %>%
  slice_min(art_pval, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(pos37, by = "rsid") %>%
  transmute(locus_label, arteriopathy, rsid, chr, pos, cad_risk_allele,
            cad_z = abs(cad_beta / cad_se), art_z = art_z_aligned, cad_p = cad_pval,
            art_p = art_pval, direction, cad_beta, art_beta_aligned, needs_manual_check) %>%
  add_hg38() %>%
  arrange(factor(arteriopathy, levels = art_order), locus_label)
write_tsv(directionality, file.path(out, "directionality_cad.tsv"))

panel_stats <- directionality %>%
  group_by(arteriopathy) %>%
  summarise(n_total = n(), n_concordant = sum(direction == "concordant"),
            n_discordant = sum(direction == "discordant"), .groups = "drop") %>%
  mutate(binom_p = map2_dbl(n_total, n_discordant, function(n, d) {
    binom.test(max(d, n - d), n, 0.5, alternative = "greater")$p.value
  }))
write_tsv(panel_stats, file.path(out, "directionality_panel_stats.tsv"))

# Figure ----
# a: loci colocalising with groups of related traits
trait_groups <- tribble(
  ~group_label, ~theme, ~pattern,
  "Coronary artery disease", "Cardiovascular",
  paste0("coronary artery disease|atherosclerotic heart|myocardial infarction|heart attack|",
         "angina|ischaemic heart|ischemic heart|family history.*heart disease|",
         "illnesses of.*heart disease|major coronary|coronary heart disease|",
         "vascular/heart problems|operation coronary|operative.*coronary|",
         "antithrombotic|salicylic|vasodilator|medication.*cholesterol.*heart|statin"),
  "Hypertension", "Cardiovascular",
  "^hypertension|high blood pressure|blood pressure medication|antihypertensive",
  "Systolic blood pressure", "Cardiovascular", "systolic blood pressure",
  "Diastolic blood pressure", "Cardiovascular", "diastolic blood pressure",
  "Pulse pressure", "Cardiovascular", "pulse pressure",
  "Lung function", "Respiratory", "fev1|fvc|peak expiratory|lung function",
  "Body fat percentage", "Anthropometric", "body fat percentage|trunk fat|arm fat|leg fat",
  "Lean mass", "Anthropometric", "lean mass|appendicular lean",
  "Height", "Anthropometric", "^standing height$|^height$|comparative height|sitting height",
  "Migraine / headache", "Neurological",
  "migraine|sumatriptan|antimigraine|anti-migraine|^headache|chronic headache|pain type.*headache"
)
theme_order <- c("Cardiovascular", "Respiratory", "Anthropometric", "Neurological")

panel_a_data <- complex %>%
  mutate(group_label = first_match(trait_name, trait_groups$pattern, trait_groups$group_label,
                                   ignore_case = TRUE)) %>%
  filter(!is.na(group_label)) %>%
  group_by(group_label, arteriopathy) %>%
  summarise(n_loci = n_distinct(coloc_group_id), .groups = "drop") %>%
  left_join(trait_groups %>% select(group_label, theme), by = "group_label") %>%
  mutate(arteriopathy = factor(arteriopathy, levels = art_order),
         theme = factor(theme, levels = theme_order))

row_order <- panel_a_data %>%
  group_by(group_label, theme) %>%
  summarise(total = sum(n_loci), .groups = "drop") %>%
  arrange(theme, desc(total)) %>%
  pull(group_label)
panel_a_data <- panel_a_data %>% mutate(group_label = factor(group_label, levels = rev(row_order)))

divider_positions <- trait_groups %>%
  transmute(group_label, theme, y_num = match(group_label, rev(row_order))) %>%
  filter(!is.na(y_num)) %>%
  arrange(y_num) %>%
  filter(!is.na(lead(theme)) & theme != lead(theme)) %>%
  pull(y_num) + 0.5

# Loci with at least one complex-trait colocalization, per arteriopathy
totals <- complex %>%
  group_by(arteriopathy) %>%
  summarise(total_loci = n_distinct(coloc_group_id), .groups = "drop")
art_labels <- setNames(sprintf("%s\n(%d loci)", totals$arteriopathy, totals$total_loci), totals$arteriopathy)

panel_a <- ggplot(panel_a_data, aes(arteriopathy, group_label)) +
  geom_hline(yintercept = divider_positions, colour = "grey88", linewidth = 0.3) +
  geom_point(aes(size = n_loci), fill = "#B2182B", shape = 21, stroke = 0.3, colour = "grey30",
             alpha = 0.85) +
  scale_size_continuous(range = c(2, 12), name = "Colocalizing\nloci", breaks = c(1, 5, 10, 20, 50)) +
  scale_x_discrete(labels = art_labels) +
  scale_y_discrete() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(face = "bold", size = 11, colour = "grey20"),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "right",
        legend.key.height = unit(0.5, "cm"), plot.margin = margin(10, 15, 10, 5))

# b: CAD-aligned z-scores at loci genome-wide significant in both traits
trait_colours <- c(AA = "#E41A1C", IA = "#377EB8", SCAD = "#4DAF4A", FMD = "#984EA3", CeAD = "#FF7F00")
label_positions <- directionality %>%
  group_by(locus_label) %>%
  summarise(label_x = mean(cad_z), label_y = mean(art_z), n_points = n(), .groups = "drop")
leader_segments <- directionality %>%
  inner_join(label_positions, by = "locus_label") %>%
  filter(n_points > 1)
x_max <- ceiling(max(directionality$cad_z) * 1.05 / 5) * 5
y_abs <- ceiling(max(abs(directionality$art_z)) * 1.05 / 5) * 5

panel_b <- ggplot(directionality, aes(cad_z, art_z)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", colour = "grey75", linewidth = 0.3,
              alpha = 0.5) +
  geom_segment(data = leader_segments, aes(x = cad_z, y = art_z, xend = label_x, yend = label_y),
               colour = "grey60", linewidth = 0.2, alpha = 0.6, inherit.aes = FALSE) +
  geom_point(aes(fill = arteriopathy), shape = 21, colour = "grey20", stroke = 0.4, size = 4.5,
             alpha = 0.9) +
  geom_text_repel(data = label_positions, aes(x = label_x, y = label_y, label = locus_label),
                  size = 4, segment.size = 0.2, segment.colour = "grey60", box.padding = 0.6,
                  max.overlaps = Inf, min.segment.length = 0.1,
                  xlim = c(-x_max * 0.05, x_max - 0.5), ylim = c(-y_abs + 0.5, y_abs - 0.5),
                  inherit.aes = FALSE) +
  scale_fill_manual(values = trait_colours, guide = "none") +
  scale_x_continuous(breaks = seq(0, x_max, 5)) +
  scale_y_continuous(breaks = seq(-y_abs, y_abs, 5)) +
  coord_cartesian(xlim = c(-1, x_max), ylim = c(-y_abs, y_abs), clip = "off") +
  labs(x = "Z-score in CAD (risk-allele aligned)", y = "Z-score in arteriopathy") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey95", linewidth = 0.2),
        legend.position = "none", plot.margin = margin(10, 10, 10, 10),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"))

# c: BP-colocalising loci shared by two or more arteriopathies
bp_groups <- c(
  SBP = "systolic blood pressure",
  DBP = "diastolic blood pressure",
  PP = "pulse pressure",
  Hypertension = paste0("^hypertension|high blood pressure|blood pressure medication|",
                        "antihypertensive|medication.*blood pressure")
)
bp_locus_data <- complex %>%
  mutate(bp_group = first_match(trait_name, bp_groups, names(bp_groups), ignore_case = TRUE)) %>%
  filter(!is.na(bp_group)) %>%
  mutate(gene_label = gene_label(display_snp)) %>%
  group_by(gene_label, arteriopathy) %>%
  summarise(n_bp_categories = n_distinct(bp_group), .groups = "drop") %>%
  group_by(gene_label) %>%
  filter(n_distinct(arteriopathy) >= 2) %>%
  mutate(n_art = n_distinct(arteriopathy), total_cat = sum(n_bp_categories)) %>%
  ungroup()
locus_order <- bp_locus_data %>%
  distinct(gene_label, n_art, total_cat) %>%
  arrange(n_art, total_cat) %>%
  pull(gene_label)
bp_locus_data <- bp_locus_data %>%
  mutate(gene_label = factor(gene_label, levels = locus_order),
         arteriopathy = factor(arteriopathy, levels = art_order))

panel_c <- ggplot(bp_locus_data, aes(arteriopathy, gene_label)) +
  geom_point(aes(size = n_bp_categories), fill = "#2166AC", shape = 21, stroke = 0.3,
             colour = "grey30", alpha = 0.9) +
  scale_size_continuous(range = c(2, 12), name = "BP categories\ncolocalizing", breaks = 1:4,
                        limits = c(1, 4)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(face = "bold", size = 11, colour = "grey20"),
        legend.title = element_text(size = 11), legend.text = element_text(size = 10),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "right", legend.key.height = unit(0.6, "cm"),
        plot.margin = margin(10, 15, 10, 5))

# d: IVW estimates for SBP and DBP
bp_colours <- c(SBP = "#B2182B", DBP = "#2166AC")
ivw <- mr_all %>%
  filter(method == "Inverse variance weighted", bp_trait %in% c("SBP", "DBP")) %>%
  mutate(or = exp(b * bp_scale),
         or_lci = exp((b - 1.96 * se) * bp_scale),
         or_uci = exp((b + 1.96 * se) * bp_scale),
         bp_trait = factor(bp_trait, levels = c("SBP", "DBP")),
         y_pos = as.numeric(factor(arteriopathy, levels = rev(art_order))) +
           if_else(bp_trait == "SBP", 0.18, -0.18),
         ci_label = sprintf("%.2f (%.2f–%.2f)", or, or_lci, or_uci))

x_breaks <- c(0.8, 1, 1.5, 2, 3, 4, 5)
x_breaks <- x_breaks[x_breaks >= min(ivw$or_lci) * 0.95 & x_breaks <= max(ivw$or_uci) * 1.05]
header_y <- length(art_order) + 0.8

panel_d <- ggplot(ivw) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = or_lci, xmax = or_uci, y = y_pos, colour = bp_trait), height = 0,
                 linewidth = 0.7) +
  geom_point(aes(x = or, y = y_pos, colour = bp_trait, fill = bp_trait), shape = 22, size = 3.2,
             stroke = 0.4) +
  geom_text(aes(x = max(x_breaks) * 1.04, y = y_pos, label = ci_label), colour = "black",
            hjust = 0, size = 3.5, show.legend = FALSE) +
  annotate("text", x = min(x_breaks) * 0.55, y = header_y, label = "Arteriopathy",
           fontface = "bold", size = 4, hjust = 1) +
  annotate("text", x = max(x_breaks) * 1.04, y = header_y, label = "OR (95% CI)",
           fontface = "bold", size = 4, hjust = 0) +
  scale_colour_manual(values = bp_colours, name = NULL) +
  scale_fill_manual(values = bp_colours, name = NULL, guide = "none") +
  scale_x_continuous(trans = "log10", breaks = x_breaks, labels = sprintf("%.1f", x_breaks),
                     expand = expansion(mult = c(0.01, 0.22))) +
  scale_y_continuous(breaks = seq_along(art_order), labels = rev(art_order),
                     expand = expansion(add = c(0.6, 0.8))) +
  coord_cartesian(clip = "off") +
  labs(x = x_lab, y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4, shape = 22, linewidth = 0,
                                                   fill = bp_colours, colour = "grey30"),
                               nrow = 1)) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(size = 11, face = "bold", colour = "black", hjust = 1),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(t = 10)),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.x = element_line(linewidth = 0.5, colour = "black"),
        axis.ticks.length.x = unit(3, "pt"), panel.grid = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 11, colour = "black"),
        legend.key = element_blank(), plot.margin = margin(15, 30, 10, 80))

row_grid <- function(a, b, labels, widths) {
  plot_grid(a, b, ncol = 2, rel_widths = widths, labels = labels, label_size = 16,
            label_fontface = "bold", label_x = 0.01, label_y = 0.99, align = "h", axis = "tb")
}
fig <- plot_grid(row_grid(panel_a, panel_b, c("a", "b"), c(1, 1)),
                 row_grid(panel_c, panel_d, c("c", "d"), c(1, 1.1)),
                 ncol = 1, rel_heights = c(1, 0.95))
ggsave(file.path(res_dir("figures"), "Figure3.pdf"), fig, width = 14, height = 11)

loci_n <- count(locus_table, arteriopathy)
write_counts("gpmap", c(
  setNames(loci_n$n, paste0("loci_", loci_n$arteriopathy)),
  setNames(totals$total_loci, paste0("loci_complex_coloc_", totals$arteriopathy)),
  complex_traits = n_distinct(complex$trait_name),
  bp_loci_2plus_arteriopathies = n_distinct(bp_locus_data$gene_label),
  bp_mr_ivw_p05 = sum(ivw$pval < 0.05),
  directionality_pairs = nrow(directionality),
  setNames(panel_stats$n_discordant, paste0("discordant_", panel_stats$arteriopathy)),
  setNames(panel_stats$n_concordant, paste0("concordant_", panel_stats$arteriopathy))
))
