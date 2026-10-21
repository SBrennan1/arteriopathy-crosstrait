# Table 1: genes associated with >= 3 arteriopathies through ASSET loci, pQTL MR or
# arterial eQTL MR, with colocalization, monogenic disease and drug annotations.

source("R/setup.R")
out <- res_dir("convergence")
p <- cfg$convergence

hits <- mr_hits()
ev <- convergence_evidence(hits$pqtl, hits$eqtl)

# Druggability ----
tract <- read_stage("drug_targets", "drug_tractability.tsv")
drug_details <- read_stage("drug_targets", "drug_details.tsv")

# Drugs whose mechanism acts on the gene product
target_drugs <- drug_details %>%
  filter(!is.na(mechanism_of_action), !is.na(drug_name),
         str_detect(tolower(mechanism_of_action), paste(p$target_engage_patterns, collapse = "|")),
         !toupper(drug_name) %in% toupper(p$drug_blacklist))

drug_status <- target_drugs %>%
  summarise(max_phase = suppressWarnings(max(clinical_phase, na.rm = TRUE)),
            approved = paste(unique(drug_name[which(clinical_phase == 4)]), collapse = "; "),
            phase3 = paste(unique(drug_name[which(clinical_phase == 3)]), collapse = "; "),
            .by = gene) %>%
  mutate(max_phase = if_else(is.infinite(max_phase), NA_real_, max_phase),
         across(c(approved, phase3), ~ na_if(.x, "")))

drug_annot <- tract %>%
  select(gene, tractability_sm, tractability_ab, tractability_protac, tractability_oc) %>%
  left_join(drug_status, by = "gene") %>%
  transmute(gene,
            tractability = case_when(tractability_sm & tractability_ab ~ "SM / AB",
                                     tractability_sm ~ "SM",
                                     tractability_ab ~ "AB",
                                     tractability_protac ~ "PROTAC",
                                     tractability_oc ~ "OC"),
            drug_status = case_when(!is.na(approved) ~ paste0("Approved (", approved, ")"),
                                    !is.na(phase3) ~ paste0("Phase 3 (", phase3, ")"),
                                    max_phase >= 1 ~ paste0("Phase ", max_phase)))

# Table 1 ----
monogenic <- read_tsv(p$monogenic, show_col_types = FALSE)
filled <- "●"   # ASSET locus or colocalized MR
open_dot <- "○" # MR only
dash <- "—"

table1 <- ev$genes %>%
  filter(in_table1) %>%
  left_join(monogenic, by = "gene") %>%
  left_join(drug_annot, by = "gene") %>%
  transmute(Gene = gene,
            GWAS = if_else(ASSET, filled, ""),
            `pQTL MR` = case_when(pQTL_coloc ~ filled, pQTL_MR ~ open_dot, TRUE ~ ""),
            `eQTL MR` = case_when(eQTL_coloc ~ filled, eQTL_MR ~ open_dot, TRUE ~ ""),
            Arteriopathies = arteriopathies,
            `Monogenic disease` = coalesce(monogenic_disease, dash),
            Tractability = coalesce(tractability, dash),
            `Drug status` = coalesce(drug_status, dash)) %>%
  arrange(desc(str_count(GWAS, filled) + str_count(`pQTL MR`, filled) + str_count(`eQTL MR`, filled)),
          desc(str_count(Arteriopathies, ",")), Gene)
write_tsv(table1, file.path(out, "table1.tsv"))

ev$genes %>%
  mutate(table1_rank = match(gene, table1$Gene)) %>%
  arrange(table1_rank, desc(n_traits), gene) %>%
  write_tsv(file.path(out, "gene_evidence.tsv"))

g <- ev$genes
write_counts("convergence", c(
  candidate_genes = nrow(g),
  table1_genes = nrow(table1),
  table1_asset = sum(table1$GWAS == filled),
  table1_pqtl_coloc = sum(table1$`pQTL MR` == filled),
  table1_pqtl_mr_only = sum(table1$`pQTL MR` == open_dot),
  table1_eqtl_coloc = sum(table1$`eQTL MR` == filled),
  table1_eqtl_mr_only = sum(table1$`eQTL MR` == open_dot),
  table1_monogenic = sum(table1$`Monogenic disease` != dash),
  table1_tractable = sum(table1$Tractability != dash),
  table1_with_drug = sum(table1$`Drug status` != dash),
  genes_2plus_modalities = sum(g$n_modalities >= 2)
))
