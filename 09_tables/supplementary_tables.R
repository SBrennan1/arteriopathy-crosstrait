# Supplementary Tables 1-31 (results/tables/supplementary_tables.xlsx). Variants are
# given by rsID with GRCh37 and GRCh38 positions. Tables 1, 11 and 12 are read from the
# files under supplementary_tables: in config.yml.

source("R/setup.R")
library(openxlsx)
out <- res_dir("tables")
fdr <- cfg$mr$fdr
strong <- cfg$qtl_coloc$strong
st <- cfg$supplementary_tables %||% list()

# Load results ----
rg <- read_result("ldsc/rg_pairwise.tsv")
asset <- read_result("asset/asset_loci.tsv")
asset_v2g <- read_result("v2g/asset_v2g.tsv")
cpassoc <- read_result("cpassoc/cpassoc_loci.tsv")
cpassoc_v2g <- read_result("v2g/cpassoc_v2g.tsv")
mtag <- read_result("mtag/mtag_leads.tsv")
mtag_v2g <- read_result("v2g/mtag_v2g.tsv")
cs <- read_result("coloc/asset_credible_sets.tsv")
coloc <- read_result("coloc/cross_trait_coloc.tsv")
pqtl <- read_result("qtl_coloc/pqtl_mr_coloc.tsv")
eqtl <- read_result("qtl_coloc/eqtl_mr_coloc.tsv")

# Downstream tables used with their own columns
other <- map(c(
  MAGMA_genes = "enrichment/magma_asset_genes.tsv",
  MAGMA_gene_sets = "enrichment/magma_asset_pathways_canonical.tsv",
  Tissue_DEG_enrichment = "enrichment/gene2func_tissue_deg.tsv",
  GO_cellular_component = "enrichment/gene2func_go_cc.tsv",
  GPMap_traits = "gpmap/gpmap_locus_annotations.tsv",
  GPMap_BP_loci = "gpmap/gpmap_bp_loci.tsv",
  BP_MR = "gpmap/bp_mr_all_methods.tsv",
  pQTL_eQTL_concordance = "qtl_coloc/pqtl_eqtl_concordance.tsv",
  MR_pathway_enrichment = "enrichment/mr_pathway_enrichment.tsv",
  AAA_target_genes = "vascular/aaa_target_genes.tsv"
), read_result)

# Tables prepared outside the pipeline (tsv, csv or xlsx)
read_external <- function(path, cols) {
  if (is.null(path) || !file.exists(path)) {
    message("not found, written as an empty table: ", path %||% "(path not set in config.yml)")
    return(as_tibble(set_names(rep(list(character()), length(cols)), cols)))
  }
  d <- switch(tolower(tools::file_ext(path)),
              xlsx = as_tibble(openxlsx::read.xlsx(path, check.names = FALSE, sep.names = " ")),
              csv = read_csv(path, show_col_types = FALSE),
              read_tsv(path, show_col_types = FALSE))
  names(d) <- str_replace_all(names(d), "_", " ")
  d
}

# Variant positions ----
pqtl_ivs <- imap_dfr(cfg$mr$pqtl_instruments, function(path, platform) {
  d <- vroom(path, show_col_types = FALSE,
             col_types = cols(rsid = "c", effect_allele = "c", other_allele = "c"))
  name_col <- intersect(c("exposure_name", "protein", "exposure"), names(d))[1]
  d %>% transmute(platform, protein = .data[[name_col]], rsid = as.character(rsid),
                  effect_allele, other_allele, eaf = suppressWarnings(as.numeric(eaf)),
                  beta = as.numeric(beta), se = as.numeric(se), pval = as.numeric(pval))
})

needs_positions <- function(d) "rsid" %in% names(d) && !all(c("pos", "pos_hg38") %in% names(d))
other_rsids <- unlist(map(keep(other, needs_positions), "rsid"))

pos_cols <- function(d) {
  transmute(d, rsID = rsid, CHR = as.character(chr), POS_GRCh37 = as.integer(pos),
            CHR_GRCh38 = as.character(chr_hg38), POS_GRCh38 = as.integer(pos_hg38))
}
positions <- list(asset, cpassoc, mtag, cs, read_result("mr/instruments_olink.tsv")) %>%
  map_dfr(pos_cols) %>%
  distinct(rsID, .keep_all = TRUE)
missing <- setdiff(na.omit(unique(c(pqtl_ivs$rsid, asset_v2g$top_snp, asset_v2g$coding_snp,
                                    coloc$hit1, coloc$hit2, other_rsids))), positions$rsID)
if (length(missing)) {
  positions <- bind_rows(positions, snp_positions(missing) %>% add_hg38() %>% pos_cols())
}

variant_id <- function(rsid, chr, pos, chr_hg38, pos_hg38) {
  tibble(rsID = rsid, CHR = as.character(chr), POS_GRCh37 = as.integer(pos),
         CHR_GRCh38 = as.character(chr_hg38), POS_GRCh38 = as.integer(pos_hg38))
}

# rsID in column `col` followed by its GRCh37 and GRCh38 positions
add_positions <- function(d, col) {
  pos <- positions %>%
    transmute(rsID, "{col}_POS_GRCh37" := POS_GRCh37, "{col}_POS_GRCh38" := POS_GRCh38)
  d %>%
    left_join(pos, by = setNames("rsID", col)) %>%
    relocate(starts_with(paste0(col, "_POS_")), .after = all_of(col))
}

# Column names in workbook style, with positions added where variants are given by rsID
format_other <- function(d) {
  if (needs_positions(d)) {
    d <- d %>%
      select(-any_of(c("chr", "pos", "chr_hg38", "pos_hg38"))) %>%
      left_join(positions %>% rename(rsid = rsID), by = "rsid")
  }
  fixed <- c(rsid = "rsID", chr = "CHR", pos = "POS_GRCh37", chr_hg38 = "CHR_GRCh38",
             pos_hg38 = "POS_GRCh38")
  nm <- names(d)
  names(d) <- ifelse(nm %in% names(fixed), fixed[nm], gsub("(^|_)([a-z])", "\\1\\U\\2", nm, perl = TRUE))
  if ("rsID" %in% names(d)) {
    d <- relocate(d, any_of(c("CHR", "POS_GRCh37", "CHR_GRCh38", "POS_GRCh38")), .after = rsID)
  }
  d
}
other <- map(other, format_other)

per_trait <- function(d, prefix) select(d, all_of(paste0(prefix, traits)))
yes_no <- function(x) if_else(x, "Yes", "No")

# Study characteristics (S1) ----
s1_cols <- c("GWAS summary statistics", "Sample Size", "Cases", "Controls", "Ancestry",
             "Case definition", "First author (year)", "PMID/DOI", "Download link")
s1 <- read_external(st$study_characteristics %||% "config/study_characteristics.tsv", s1_cols)

# LDSC (S2, S3) ----
s2 <- read_result("ldsc/h2_liability.tsv") %>%
  transmute(Trait = trait, N_SNPs = n_snps,
            h2_Observed = round(h2_obs, 3), h2_Observed_SE = round(h2_obs_se, 3),
            h2_Liability = round(h2_liab, 3), h2_Liability_SE = round(h2_liab_se, 3),
            Sample_Prevalence = round(samp_prev, 3), Population_Prevalence = pop_prev,
            Intercept = round(intercept, 3), Intercept_SE = round(intercept_se, 3))

# BH across the ten trait pairs, as in the heatmap
s3 <- rg %>%
  transmute(`Trait 1` = std_trait(p1), `Trait 2` = std_trait(p2), rg = round(rg, 2), SE = round(se, 2),
            Z = round(z, 2), P = signif(p, 2), `FDR-adjusted P` = signif(p.adjust(p, method = "BH"), 2),
            `Genetic covariance intercept` = round(gcov_int, 3), `Intercept SE` = round(gcov_int_se, 3))

# Intercept sensitivity (S4, S6)
intercept_check <- function(file) {
  read_result(file.path("sensitivity", file)) %>%
    transmute(rsid = locus_rsid, P_Intercept_Corrected = p_rescaled,
              Traits_Intercept_Corrected = traits_rescaled,
              Retained_Intercept_Corrected = yes_no(retained))
}

# ASSET loci and genes (S4, S5) ----
s4 <- asset %>%
  left_join(intercept_check("intercept_asset.tsv"), by = "rsid") %>%
  left_join(read_result("asset/asset_cytobands.tsv") %>% select(rsid, cytoband), by = "rsid") %>%
  transmute(variant_id(rsid, chr, pos, chr_hg38, pos_hg38), Cytoband = cytoband, A1, A2,
            Nearest_Gene = nearest_gene, ASSET_P = asset_p,
            Subset_Positive = asset_subset_pos, Subset_Negative = asset_subset_neg,
            N_Traits_Nominal = n_nominal, Traits_Positive = traits_pos,
            Traits_Negative = traits_neg, Effect_Pattern = effect_pattern,
            P_Intercept_Corrected, Traits_Intercept_Corrected, Retained_Intercept_Corrected) %>%
  bind_cols(per_trait(asset, "P_"), per_trait(asset, "Z_")) %>%
  arrange(ASSET_P)

s5 <- asset_v2g %>%
  transmute(variant_id(locus_rsid, chr, pos, chr_hg38, pos_hg38),
            Gene = gene, Mapping_Evidence = evidence,
            Top_SNP = top_snp, Top_PIP = top_pip, N_Credible_Set_SNPs = n_cs_snps,
            Nearest_Gene = nearest_gene, Nearest_Distance = nearest_distance,
            Coding_SNP = coding_snp, Coding_Gene = coding_gene, Coding_PIP = coding_pip,
            Coding_Consequence = coding_consequence,
            Coloc_Gene = coloc_gene, Coloc_Tissue = coloc_tissue, Coloc_PP_H4 = coloc_pp_h4,
            Coloc_Trait = coloc_trait,
            Tiered_Gene = gene_tiered, Tiered_Evidence = evidence_tiered, Notes = override_note) %>%
  add_positions("Top_SNP") %>%
  add_positions("Coding_SNP")

# CPASSOC cross-trait loci (S6) ----
s6 <- cpassoc %>%
  left_join(intercept_check("intercept_cpassoc.tsv"), by = "rsid") %>%
  left_join(cpassoc_v2g %>%
              select(rsid, gene, evidence, vep_consequence, coloc_gene, coloc_tissue,
                     coloc_pp_h4, coloc_trait, override_note),
            by = "rsid") %>%
  transmute(variant_id(rsid, chr, pos, chr_hg38, pos_hg38), A1, A2,
            P_SHet, P_SHom, Cross_Trait_Tier = cross_trait_tier,
            Gene = gene, Mapping_Evidence = evidence, In_ASSET_Genes = yes_no(gene %in% asset_v2g$gene),
            VEP_Consequence = vep_consequence,
            Nearest_Gene = nearest_gene, Nearest_Distance = nearest_distance,
            Coloc_Gene = coloc_gene, Coloc_Tissue = coloc_tissue, Coloc_PP_H4 = coloc_pp_h4,
            Coloc_Trait = coloc_trait, Notes = override_note,
            N_Traits_Nominal = n_nominal, Traits_Positive = traits_pos,
            Traits_Negative = traits_neg, Effect_Pattern = effect_pattern,
            P_Intercept_Corrected, Traits_Intercept_Corrected, Retained_Intercept_Corrected) %>%
  bind_cols(per_trait(cpassoc, "P_"), per_trait(cpassoc, "Z_")) %>%
  arrange(P_SHet)

# MTAG (S7, S8, S9) ----
s7 <- read_result("mtag/mtag_parameters.tsv") %>%
  transmute(MTAG_Pair = mtag_pair, Trait = trait, N = n, rg_MTAG = round(rg_mtag, 3),
            GWAS_chi2 = round(gwas_chi2, 3), MTAG_chi2 = round(mtag_chi2, 3),
            chi2_change = round(chi2_change, 3), GWAS_equiv_N = gwas_equiv_n,
            maxFDR = round(max_fdr, 3))

s8 <- mtag %>%
  left_join(mtag_v2g %>% select(rsid, gene, evidence, override_note), by = "rsid") %>%
  transmute(variant_id(rsid, chr, pos, chr_hg38, pos_hg38),
            Gene = gene, Mapping_Evidence = evidence, Trait = target_trait, MTAG_Pair = mtag_pair,
            Effect_Allele = effect_allele_univ, Other_Allele = other_allele_univ,
            Beta_Univariate = round(beta_univ, 4), P_Univariate = signif(p_univ, 3),
            Beta_MTAG = round(beta_mtag, 4), P_MTAG = signif(p_mtag, 3),
            Uplift = round(uplift, 2), Novel_MTAG = yes_no(novel_mtag),
            Previously_Reported = yes_no(previously_reported), Truly_Novel = yes_no(novel),
            Mapping_Notes = override_note) %>%
  arrange(Trait, P_MTAG)

s9 <- s8 %>%
  filter(Truly_Novel == "Yes") %>%
  select(-Novel_MTAG, -Previously_Reported, -Truly_Novel, -Mapping_Notes)

# Cross-trait colocalization (S10) ----
s10 <- coloc %>%
  left_join(asset_v2g %>% select(locus_rsid, gene), by = "locus_rsid") %>%
  transmute(variant_id(locus_rsid, chr, pos, chr_hg38, pos_hg38),
            Gene = gene, Nearest_Gene = nearest_gene, Trait_1 = trait1, Trait_2 = trait2,
            N_SNPs = n_snps,
            PP_H0 = PP.H0.abf, PP_H1 = PP.H1.abf, PP_H2 = PP.H2.abf, PP_H3 = PP.H3.abf,
            PP_H4 = PP.H4.abf, Class = class,
            SuSiE_PP_H3 = susie_PP.H3, SuSiE_PP_H4 = susie_PP.H4, SuSiE_Class = susie_class,
            Hit_1 = hit1, Hit_2 = hit2) %>%
  add_positions("Hit_1") %>%
  add_positions("Hit_2")

# Cross-ancestry evaluation (S11, S12) ----
s11_cols <- c("Mapped Gene", "rsID", "CHR", "POS (GRCh37)", "POS (GRCh38)", "Effect Allele",
              "BETA EUR", "P EUR", "BETA MTAG", "P MTAG", "BETA EAS", "P EAS")
s12_cols <- c("Gene", "Trait", "rsID", "External ancestry", "Effect Allele", "EAF EUR", "BETA EUR",
              "EUR P", "External EAF", "BETA External", "External P", "Concordant", "P<0.0083")
s11 <- read_external(st$cross_ancestry_mtag, s11_cols)
s12 <- read_external(st$cross_ancestry_coloc, s12_cols)

# Negative control phenotypes (S13) ----
s13 <- read_tsv(cfg$negative_controls$draws, show_col_types = FALSE) %>%
  summarise(`Matched arteriopathy` = paste(sort(unique(trait)), collapse = ", "),
            `Genetic correlation` = if_else(any(ldsc_subset), "Yes", "No"),
            .by = c(phenocode, description, category, n_cases_EUR, n_controls_EUR)) %>%
  transmute(Phecode = phenocode, Description = description, Category = category,
            `Cases (EUR)` = n_cases_EUR, `Controls (EUR)` = n_controls_EUR,
            `Matched arteriopathy`, `Genetic correlation`) %>%
  arrange(`Matched arteriopathy`, Description)

# MR (S21, S22, S23, S24) ----
s21 <- pqtl_ivs %>%
  left_join(positions, by = c(rsid = "rsID")) %>%
  transmute(Platform = recode(platform, Olink = "Olink (UKB-PPP)", deCODE = "SomaScan (deCODE)"),
            Protein = protein, rsID = rsid, CHR, POS_GRCh37, CHR_GRCh38, POS_GRCh38,
            Effect_allele = effect_allele, Other_allele = other_allele, EAF = eaf,
            Beta = beta, SE = se, P = pval, F_statistic = (beta / se)^2) %>%
  arrange(Platform, Protein, P)

or_cols <- function(d) {
  mutate(d, OR = exp(Beta), OR_95CI_lower = exp(Beta - 1.96 * SE),
         OR_95CI_upper = exp(Beta + 1.96 * SE), .after = SE)
}
round_mr <- function(d) {
  mutate(d, across(c(Beta, SE, OR, OR_95CI_lower, OR_95CI_upper), ~ round(.x, 3)),
         across(starts_with("PP_H4"), ~ round(.x, 2)),
         P = signif(P, 3), FDR_q = signif(FDR_q, 3))
}

# Coloc_Validated: PP.H4 >= 0.8 (pQTL by ABF or SuSiE)
s22 <- pqtl %>%
  filter(qval < fdr) %>%
  transmute(Gene = exposure, Trait = outcome, Platform = platform, Beta = b, SE = se,
            P = pval, FDR_q = qval, N_Instruments = nsnp, Instruments = snps,
            PP_H4_ABF = abf_pph4, PP_H4_SuSiE = susie_pph4,
            Coloc_Validated = case_when(abf_pph4 >= strong | susie_pph4 >= strong ~ "Yes",
                                        !is.na(abf_pph4) | !is.na(susie_pph4) ~ "No",
                                        TRUE ~ "Not tested"),
            Coloc_Evidence = coloc_evidence) %>%
  or_cols() %>%
  round_mr() %>%
  arrange(FDR_q, P)

s23 <- read_result("mr/pqtl_replication.tsv") %>%
  left_join(pqtl %>% select(platform, protein = exposure, outcome, coloc_evidence),
            by = c("platform", "protein", "outcome")) %>%
  transmute(Protein = protein, Trait = outcome, Discovery_Platform = platform,
            Discovery_OR = exp(disc_b), Discovery_FDR_q = disc_q, Coloc_Evidence = coloc_evidence,
            Replication_Study = recode(study, Fenland = "Fenland (SomaScan)", SCALLOP = "SCALLOP (Olink)"),
            Method = rep_method, N_Instruments = rep_nsnp, Instruments = rep_snps,
            Replication_OR = exp(rep_b), Replication_OR_95CI_lower = exp(rep_b - 1.96 * rep_se),
            Replication_OR_95CI_upper = exp(rep_b + 1.96 * rep_se), Replication_P = rep_p,
            Same_Direction = if_else(tested, yes_no(concordant), NA_character_),
            Replicated = if_else(tested, yes_no(replicated), paste0("Not tested (", tolower(not_tested), ")")))

s24 <- eqtl %>%
  filter(qval < fdr) %>%
  transmute(Gene = exposure, Gene_ID = gene_id, Trait = outcome, Tissue = gsub("_", " ", tissue),
            Beta = b, SE = se, P = pval, FDR_q = qval, N_Instruments = nsnp, Instruments = snps,
            PP_H4 = PP.H4,
            Coloc_Validated = case_when(PP.H4 >= strong ~ "Yes", !is.na(PP.H4) ~ "No",
                                        TRUE ~ "Not tested"),
            Coloc_Evidence = coloc_evidence) %>%
  or_cols() %>%
  round_mr() %>%
  arrange(FDR_q, P)

# Genes associated with >= 3 arteriopathies (S27) ----
# Direction of effect across the two ASSET subsets
asset_direction <- asset_v2g %>%
  filter(!is.na(gene)) %>%
  distinct(gene, locus_rsid) %>%
  inner_join(asset %>% select(locus_rsid = rsid, asset_subset_pos, asset_subset_neg), by = "locus_rsid") %>%
  mutate(across(c(asset_subset_pos, asset_subset_neg), ~ str_replace_all(coalesce(.x, ""), ",", ", ")),
         dir = if_else(asset_subset_pos == "" | asset_subset_neg == "", "Same",
                       paste0("Opposite (", asset_subset_pos, " vs ", asset_subset_neg, ")"))) %>%
  summarise(ASSET_Direction = paste(unique(dir), collapse = "; "), .by = gene)

s27 <- read_result("convergence/gene_evidence.tsv") %>%
  filter(in_table1) %>%
  left_join(asset_direction, by = "gene") %>%
  transmute(Gene = gene, ASSET = yes_no(ASSET), ASSET_Direction, pQTL_MR = yes_no(pQTL_MR),
            pQTL_Coloc = yes_no(pQTL_coloc), eQTL_MR = yes_no(eQTL_MR),
            eQTL_Coloc = yes_no(eQTL_coloc), N_Modalities = n_modalities, N_Traits = n_traits,
            Arteriopathies = arteriopathies) %>%
  arrange(desc(N_Traits), desc(N_Modalities), Gene)

# Cell-type specificity (S28): all atlases and analyses in one table ----
enrich_rows <- function(path, atlas, analysis) {
  d <- read_result(path)
  if (!"condition" %in% names(d)) d$condition <- NA_character_
  d %>% transmute(Atlas = atlas, Analysis = analysis,
                  Donor_group = str_to_sentence(condition), Cell_type = cell_type,
                  N_target_genes = n_units, Observed_specificity = observed, Expected_specificity = null_mean,
                  Fold = observed / null_mean,
                  `Comparison_fold_95%_range` = sprintf("%.2f-%.2f", null_lo / null_mean, null_hi / null_mean),
                  Direction = if_else(Fold >= 1, "Higher", "Lower"),
                  P = if_else(Fold >= 1, p_greater, p_less),
                  P_FDR = if_else(Fold >= 1, q_greater, q_less))
}
s28 <- bind_rows(
  enrich_rows("vascular/hca_celltype_enrichment.tsv", "Arterial atlas", "Primary"),
  enrich_rows("vascular/hca_celltype_enrichment_exprmatched.tsv", "Arterial atlas", "Expression-matched comparison genes"),
  enrich_rows("vascular/hca_celltype_enrichment_nearest_gene.tsv", "Arterial atlas", "Nearest-gene assignment"),
  enrich_rows("vascular/brain_ia_specificity.tsv", "Brain aneurysm atlas", "Primary"),
  enrich_rows("vascular/brain_ia_specificity_nearest_gene.tsv", "Brain aneurysm atlas", "Nearest-gene assignment")
)

# Drug targets (S30, S31) ----
salt_words <- c("ANHYDROUS", "SODIUM", "HYCLATE", "ESYLATE", "PHOSPHATE", "MESYLATE", "HYDROCHLORIDE",
                "DIHYDROCHLORIDE", "CALCIUM", "POTASSIUM", "ACETATE", "MALEATE", "TARTRATE", "SULFATE",
                "CITRATE", "SUCCINATE", "BESYLATE", "TOSYLATE", "FUMARATE", "HYDRATE", "MONOHYDRATE")
# Drug names without salt forms; code names stay upper case
tidy_drugs <- function(x) {
  map_chr(x, function(v) {
    if (is.na(v)) return(NA_character_)
    more <- str_extract(v, "\\+\\d+\\s*$")
    parts <- str_split_1(str_remove(v, "\\s*\\+\\d+\\s*$"), ",\\s*")
    d <- str_squish(str_remove_all(toupper(parts), paste0("\\b(", paste(salt_words, collapse = "|"), ")\\b")))
    d <- unique(if_else(str_detect(d, "\\d"), d, tolower(d)))
    paste(c(paste(d, collapse = ", "), if (!is.na(more)) str_squish(more)), collapse = " ")
  })
}
phase_text <- function(x) case_when(x == 4 ~ "Approved", x >= 1 ~ paste("Phase", x))

tiers_raw <- read_result("prioritization/target_tiers.tsv")
s30 <- tiers_raw %>%
  transmute(Gene = gene, Tier = tier, Tier_1_locus = tier1_locus, Gene_type = gene_type,
            Arteriopathies_MR = mr_traits, Colocalized = coloc_traits, QTL = qtl,
            MR_estimates = mr_effects, Implied_action = implied_action,
            pQTL_replication = pqtl_replication, pQTL_protein_altering_variant = pqtl_protein_altering,
            Mendelian_disease = monogenic_disease, Mendelian_direction = mendelian_direction,
            Tractable_modalities = tractable, Compatible_modalities = feasible_modalities,
            LOEUF = round(loeuf, 2), LoF_intolerant = if_else(is.na(lof_intolerant), NA_character_, yes_no(lof_intolerant)),
            Drugs_same_direction = tidy_drugs(matching_drugs), Same_direction_max_phase = phase_text(matching_max_phase),
            Drugs_opposite_direction = tidy_drugs(opposite_drugs), Opposite_direction_max_phase = phase_text(opposite_max_phase),
            Artery_TPM = round(artery_tpm, 1), Highest_tissue = top_tissue, HCA_cell_type = hca_cell_type,
            Ensembl_ID = ensembl_id) %>%
  arrange(Tier, Gene)

# Drugs for the tiered genes, with direction relative to the implied action
action_pat <- function(x) paste(cfg$prioritization$action_patterns[[x]], collapse = "|")
s31 <- read_result("drug_targets/drug_details.tsv") %>%
  filter(!is.na(drug_name)) %>%
  inner_join(tiers_raw %>% select(gene, tier, implied_action), by = "gene") %>%
  mutate(a = tolower(coalesce(action_type, "")),
         drug_action = case_when(str_detect(a, action_pat("inhibit")) ~ "Inhibit",
                                 str_detect(a, action_pat("activate")) ~ "Activate"),
         excluded = toupper(drug_name) %in% toupper(cfg$convergence$drug_blacklist),
         Direction = case_when(excluded ~ "Excluded (does not act on the gene product)",
                               is.na(drug_action) | is.na(implied_action) | implied_action == "Opposing" ~ "Not classified",
                               drug_action == implied_action ~ "Same as implied",
                               TRUE ~ "Opposite to implied")) %>%
  transmute(Gene = gene, Tier = tier, Implied_action = implied_action, Drug = tidy_drugs(drug_name),
            ChEMBL_ID = drug_id, Drug_type = drug_type, Mechanism_of_action = mechanism_of_action,
            Action_type = str_to_sentence(action_type), Highest_phase = phase_text(clinical_phase),
            Direction, Indications = indications, phase_sort = clinical_phase) %>%
  distinct() %>%
  arrange(Tier, Gene, desc(phase_sort), Drug) %>%
  select(-phase_sort)

# Tables, titles and columns ----
# cols: the columns kept, in order, as named after formatting (underscores as spaces);
# NULL keeps every column of a table built above.
spec <- list(
  list(s1,  "Basic characteristics of included studies", s1_cols),
  list(s2,  "SNP-based heritability estimates and LDSC summary metrics for the five arteriopathies", NULL),
  list(s3,  "Pairwise genetic correlation estimates between arteriopathies from LD score regression", NULL),
  list(s4,  "ASSET cross-trait association results across non-syndromic arteriopathies", NULL),
  list(s5,  "Candidate genes for the ASSET loci", NULL),
  list(s6,  "CPASSOC cross-trait association results across non-syndromic arteriopathies", NULL),
  list(s7,  "Pairwise MTAG analysis parameters and false discovery rate estimates", NULL),
  list(s8,  "Prioritized lead variants from pairwise MTAG analyses", NULL),
  list(s9,  "Novel loci identified by MTAG", NULL),
  list(s10, "Pairwise cross-trait colocalization at ASSET loci", NULL),
  list(s11, "Cross-ancestry evaluation of MTAG-prioritized IA–AAD associations", s11_cols),
  list(s12, "Cross-ancestry evaluation of colocalized IA–AAD associations", s12_cols),
  list(s13, "Negative control phenotypes (Pan-UK Biobank) matched to each arteriopathy", NULL),
  list(other$MAGMA_genes, "MAGMA gene-level association results for cross-trait ASSET analysis",
       c("Symbol", "Ensembl Id", "CHR", "Start", "Stop", "Nsnps", "Nparam", "Zstat", "P",
         "P Bonferroni", "P Adjust Method")),
  list(other$MAGMA_gene_sets, "MAGMA gene-set enrichment results for ASSET loci",
       c("Pathway Id", "Pathway", "Database", "N Genes", "Beta", "Beta Std", "Se", "P", "P Fdr",
         "P Bonferroni")),
  list(other$Tissue_DEG_enrichment, "Tissue-specific expression enrichment of cross-trait prioritized genes",
       c("Tissue", "N Genes In Set", "N Overlapping", "P", "Adj P")),
  list(other$GO_cellular_component, "Gene Ontology cellular component enrichment of cross-trait prioritized genes",
       c("Gene Set", "N Genes In Set", "N Overlapping", "P", "Adj P", "Genes")),
  list(other$GPMap_traits, "GPMap colocalization results, molecular QTL annotations and variant-to-gene assignments for arteriopathy-associated loci",
       c("Arteriopathy", "Candidate Variant", "CHR", "POS GRCh37", "Mapped Gene", "Mapping Source",
         "Gpmap All Genes", "Gpmap Tissues", "Gpmap Data Types", "N Molecular Colocs", "N Complex Traits",
         "Trait Categories", "Top Colocalizing Traits", "Nearest Gene", "Nearest Distance",
         "Asset V2g Gene", "Asset V2g Evidence")),
  list(other$GPMap_BP_loci, "Loci showing colocalizing associations with blood pressure-related phenotypes across arteriopathies",
       c("Gene", "Candidate Variants", "N Arteriopathies", "Arteriopathies", "Bp Traits")),
  list(other$BP_MR, "Mendelian randomization results for the effect of blood pressure traits on arteriopathies",
       c("Bp Trait", "Arteriopathy", "Method", "Nsnp", "B", "Se", "Exposure Unit", "Or", "Or Lci",
         "Or Uci", "Pval")),
  list(s21, "Cis-pQTL instruments for Olink- and SomaScan-measured protein abundance in UK Biobank and deCODE", NULL),
  list(s22, "FDR-significant proteome-wide MR with colocalization", NULL),
  list(s23, "Replication of the FDR-significant pQTL MR results in Fenland (SomaScan) and SCALLOP (Olink)", NULL),
  list(s24, "FDR-significant transcriptome-wide MR with colocalization", NULL),
  list(other$pQTL_eQTL_concordance, "Direction of pQTL and eQTL MR estimates for gene-trait pairs significant in both",
       c("Platform", "Gene", "Outcome", "B Pqtl", "Pqtl Coloc", "Tissue", "B Eqtl", "Eqtl Coloc",
         "Arterial", "Concordant")),
  list(other$MR_pathway_enrichment, "Pathway enrichment of pQTL and arterial eQTL MR genes",
       c("Source", "Database", "Id", "Description", "Gene Ratio", "Bg Ratio", "Pvalue", "P Adjust",
         "Qvalue", "Genes", "Count", "Theme")),
  list(s27, "Genes associated with at least three arteriopathies through ASSET loci, pQTL MR or arterial eQTL MR", NULL),
  list(s28, "Cell-type specificity of target genes relative to other genes at the same loci", NULL),
  list(other$AAA_target_genes, "Differential expression of AAD target genes in abdominal aortic aneurysm vs donor aorta",
       c("Gene", "Sources", "Log Fc", "T", "P", "Fdr", "Ave Expr")),
  list(s30, "Genetic support, implied drug action, tractability, constraint, drugs and expression for each gene", NULL),
  list(s31, "Drugs targeting the prioritized genes (Open Targets Platform)", NULL)
)

# AA is reported as AAD (aortic aneurysm and dissection)
aad <- function(x) str_replace_all(x, "(?<![A-Za-z0-9])AA(?![A-Za-z0-9])", "AAD")
# P values to 3 and other numbers to 4 significant figures
is_p <- function(n) str_detect(n, "(^|_| )(P|p|q|Q|FDR|pval|qval)($|_| )|^P[_ ]|P$|FDR|Fdr|q_|qval|pval|Pval|Qvalue")
format_table <- function(d) {
  num <- names(d)[map_lgl(d, ~ is.numeric(.x) && any(.x != round(.x), na.rm = TRUE))]
  d <- mutate(d, across(all_of(num[is_p(num)]), ~ signif(.x, 3)),
              across(all_of(num[!is_p(num)]), ~ signif(.x, 4)),
              across(where(is.logical), ~ if_else(.x, "Yes", "No")),
              across(where(is.character), aad))
  names(d) <- aad(str_replace_all(names(d), "_", " "))
  d
}
keep_cols <- function(d, cols, title) {
  if (is.null(cols)) return(d)
  absent <- setdiff(cols, names(d))
  if (length(absent)) {
    stop("'", title, "': columns not found: ", paste(absent, collapse = ", "),
         "\n  available: ", paste(names(d), collapse = ", "), call. = FALSE)
  }
  select(d, all_of(cols))
}

# Display names for column headers
pretty <- c(
  "Ensembl Id" = "Ensembl ID", "Nsnps" = "N SNPs", "Nsnp" = "N SNPs", "Nparam" = "N parameters",
  "Zstat" = "Z", "P Adjust Method" = "P adjustment", "Pathway Id" = "Pathway ID",
  "Beta Std" = "Beta (standardized)", "Se" = "SE", "P Fdr" = "P FDR", "Adj P" = "Adjusted P",
  "Gpmap All Genes" = "GPMap genes", "Gpmap Tissues" = "GPMap tissues",
  "Gpmap Data Types" = "GPMap data types", "Bp Traits" = "BP traits", "Bp Trait" = "BP trait",
  "B" = "Beta", "Or" = "OR", "Or Lci" = "OR 95CI lower", "Or Uci" = "OR 95CI upper", "Pval" = "P",
  "B Pqtl" = "Beta pQTL", "Pqtl Coloc" = "pQTL coloc", "B Eqtl" = "Beta eQTL", "Eqtl Coloc" = "eQTL coloc",
  "Id" = "ID", "Pvalue" = "P", "P Adjust" = "P FDR", "Qvalue" = "q", "Log Fc" = "Log fold change",
  "T" = "t", "Fdr" = "P FDR", "Ave Expr" = "Average expression",
  "Asset V2g Gene" = "Assigned gene", "Asset V2g Evidence" = "Gene assignment evidence"
)
rename_pretty <- function(d) {
  i <- names(d) %in% names(pretty)
  names(d)[i] <- pretty[names(d)[i]]
  d
}

titles <- map_chr(spec, 2)
built <- map(spec, function(s) {
  if (is.null(s[[1]])) return(simpleError(paste0("table missing for '", s[[2]], "'")))
  tryCatch(s[[1]] %>% format_table() %>% keep_cols(s[[3]], s[[2]]) %>% rename_pretty(),
           error = function(e) e)
})
failed <- map_lgl(built, ~ inherits(.x, "error"))
if (any(failed)) {
  stop(paste0("Supplementary Table ", which(failed), ": ", map_chr(built[failed], conditionMessage),
              collapse = "\n\n"), call. = FALSE)
}
tables <- built
numbers <- paste0("Supplementary Table ", seq_along(tables))
sheet_names <- paste0("S", seq_along(tables))

# Write workbook ----
wb <- createWorkbook()
bold <- createStyle(textDecoration = "bold")
addWorksheet(wb, "Contents")
writeData(wb, "Contents", tibble(Table = numbers, Title = aad(titles), Sheet = sheet_names),
          headerStyle = bold)
setColWidths(wb, "Contents", cols = 1:3, widths = c(24, 120, 8))
for (i in seq_along(tables)) {
  sh <- sheet_names[i]
  addWorksheet(wb, sh)
  writeData(wb, sh, paste0(numbers[i], ". ", aad(titles[i])), startRow = 1)
  addStyle(wb, sh, bold, rows = 1, cols = 1)
  writeData(wb, sh, tables[[i]], startRow = 3, headerStyle = bold)
  freezePane(wb, sh, firstActiveRow = 4)
  setColWidths(wb, sh, cols = seq_len(max(1, ncol(tables[[i]]))), widths = "auto")
}
saveWorkbook(wb, file.path(out, "supplementary_tables.xlsx"), overwrite = TRUE)
message("Wrote ", length(tables), " supplementary tables to ", file.path(out, "supplementary_tables.xlsx"))
