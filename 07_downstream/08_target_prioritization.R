# Tiered drug target table:
#   Tier 1  MR (q < 0.05) with colocalization (PP.H4 >= 0.8) for the same exposure and
#           arteriopathy
#   Tier 2  MR without colocalization
#   Tier 3  ASSET locus gene without MR support
# Genes are annotated with the drug action implied by the MR direction, tractable
# modalities, drugs, constraint, monogenic disease, GTEx v10 expression and HCA cell type.
# Run after 07_downstream/04 and 07 and 11_vascular_expression/03.

source("R/setup.R")
out <- res_dir("prioritization")
p <- cfg$prioritization
fdr <- cfg$mr$fdr
pph4_cut <- cfg$convergence$coloc_pph4

ev <- convergence_evidence(mr_hits()$pqtl, mr_hits()$eqtl)$genes %>%
  filter(!str_detect(gene, "^ENSG"))

# MR evidence: one row per gene, arteriopathy and exposure (platform or tissue) ----
mr <- bind_rows(
  read_stage("qtl_coloc", "pqtl_mr_coloc.tsv") %>%
    filter(qval < fdr) %>%
    transmute(gene = exposure, trait = std_trait(outcome), source = "pQTL", exposure_set = platform,
              b, se, qval, pph4 = pmax(abf_pph4, susie_pph4, na.rm = TRUE)),
  read_stage("qtl_coloc", "eqtl_mr_coloc.tsv") %>%
    filter(qval < fdr, tissue %in% cfg$mr$eqtl$arterial_tissues) %>%
    transmute(gene = exposure, trait = std_trait(outcome), source = "eQTL", exposure_set = tissue,
              b, se, qval, pph4 = PP.H4)
) %>%
  filter(gene %in% ev$gene) %>%
  mutate(coloc = coalesce(pph4 >= pph4_cut, FALSE))

# One lead estimate per gene and arteriopathy: colocalized first, then pQTL, then lowest q
lead <- mr %>%
  arrange(desc(coloc), source != "pQTL", qval) %>%
  slice_head(n = 1, by = c(gene, trait))

# Implied action: inhibit if higher levels raise risk for every arteriopathy, activate if
# they lower it, otherwise opposing
gene_mr <- lead %>%
  arrange(gene, match(trait, traits)) %>%
  summarise(tier = if_else(any(coloc), 1L, 2L),
            mr_traits = paste(trait, collapse = ", "),
            coloc_traits = na_if(paste(trait[coloc], collapse = ", "), ""),
            risk_up = na_if(paste(trait[b > 0], collapse = ", "), ""),
            risk_down = na_if(paste(trait[b < 0], collapse = ", "), ""),
            implied_action = case_when(all(b > 0) ~ "Inhibit", all(b < 0) ~ "Activate", TRUE ~ "Opposing"),
            or_text = paste0(trait, " ", sprintf("%.2f", exp(b)), collapse = "; "),
            qtl = na_if(paste(sort(unique(source[coloc])), collapse = ", "), ""),
            mr_effects = paste0(trait, " OR ", sprintf("%.2f", exp(b)), " (", source, " ", exposure_set,
                                if_else(coloc, sprintf(", PP.H4 %.2f", pph4), ""), ")", collapse = "; "),
            .by = gene)

# Optional direction overrides (config: direction_overrides)
ov <- p$direction_overrides %||% list()
if (length(ov)) {
  ov <- tibble(gene = names(ov), override = map_chr(ov, "action"), direction_note = map_chr(ov, "reason"))
  gene_mr <- gene_mr %>%
    left_join(ov, by = "gene") %>%
    mutate(implied_action = coalesce(override, implied_action)) %>%
    select(-override)
} else gene_mr$direction_note <- NA_character_

# Arteriopathies where the significant exposures (platforms, tissues) disagree in sign
disagree <- mr %>%
  summarise(differs = n_distinct(sign(b)) > 1, .by = c(gene, trait)) %>%
  summarise(exposures_disagree = na_if(paste(trait[differs], collapse = ", "), ""), .by = gene)

# Drugs, by the change in target activity their mechanism causes ----
action_of <- function(a) {
  a <- tolower(coalesce(a, ""))
  pat <- function(x) paste(p$action_patterns[[x]], collapse = "|")
  case_when(str_detect(a, pat("inhibit")) ~ "Inhibit", str_detect(a, pat("activate")) ~ "Activate")
}
top_drugs <- function(name, phase, n = 3) {
  u <- unique(name[order(-phase)])
  if (length(u) == 0) return(NA_character_)
  paste0(paste(head(u, n), collapse = ", "), if (length(u) > n) paste0(" +", length(u) - n) else "")
}
max_or_na <- function(x) if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

drugs <- read_stage("drug_targets", "drug_details.tsv") %>%
  filter(!is.na(drug_name), !toupper(drug_name) %in% toupper(cfg$convergence$drug_blacklist)) %>%
  mutate(drug_action = action_of(action_type))
print(count(drugs, action_type, drug_action, sort = TRUE), n = Inf)

drug_status <- drugs %>%
  summarise(max_clinical_phase = max_or_na(clinical_phase),
            approved_drugs = top_drugs(drug_name[which(clinical_phase == 4)], clinical_phase[which(clinical_phase == 4)], n = 5),
            .by = gene)

drug_match <- drugs %>%
  inner_join(gene_mr %>% select(gene, implied_action), by = "gene") %>%
  filter(implied_action != "Opposing", !is.na(drug_action)) %>%
  mutate(matches = drug_action == implied_action) %>%
  summarise(matching_drugs = top_drugs(drug_name[matches], clinical_phase[matches]),
            matching_max_phase = max_or_na(clinical_phase[matches]),
            opposite_drugs = top_drugs(drug_name[!matches], clinical_phase[!matches]),
            opposite_max_phase = max_or_na(clinical_phase[!matches]),
            .by = gene)

# Tractability, constraint, burden and safety (Open Targets) ----
tract <- read_stage("drug_targets", "drug_tractability.tsv") %>%
  mutate(tractable = pmap_chr(list(tractability_sm, tractability_ab, tractability_protac, tractability_oc),
                              function(...) na_if(paste(c("SM", "AB", "PROTAC", "OC")[c(...)], collapse = ", "), "")))

# Monogenic disease; lof_mechanism marks loss-of-function disorders (lower level, higher risk)
mono <- read_tsv(cfg$convergence$monogenic, show_col_types = FALSE) %>%
  mutate(lof_mechanism = coalesce(as.logical(lof_mechanism), FALSE))

# GENCODE gene type and position
gencode <- fread(cmd = paste("gzip -dc", shQuote(cfg$gtex_v10$gencode), "| awk -F'\\t' '$3==\"gene\"'"),
                 header = FALSE, sep = "\t", select = c(1, 4, 9), col.names = c("chr", "start", "attr"))[!grepl("_PAR_Y", attr)]
gtypes <- gencode[, .(ensembl_id = sub('.*gene_id "(ENSG\\d+).*', "\\1", attr),
                      chr = sub("^chr", "", chr), start,
                      gene = sub('.*gene_name "([^"]+)".*', "\\1", attr),
                      gene_type = sub('.*gene_type "([^"]+)".*', "\\1", attr))]

# pQTL instruments tagging a protein-altering variant in their own protein
pav_file <- file.path(cfg$results_dir, "prioritization/pqtl_protein_altering.tsv")
pav <- if (file.exists(pav_file)) {
  read_tsv(pav_file, show_col_types = FALSE) %>%
    filter(protein_altering) %>%
    summarise(pqtl_protein_altering = paste(unique(sprintf("%s (%s, r2 %.2f)", pav_variant, pav_change, pav_r2)),
                                            collapse = "; "), .by = gene)
} else {
  message("no protein-altering variant check (run 07_downstream/07_pqtl_protein_altering.R first)")
  tibble(gene = character(), pqtl_protein_altering = character())
}

# pQTL replication
repl_file <- file.path(cfg$results_dir, "mr/pqtl_replication.tsv")
repl <- if (file.exists(repl_file)) {
  read_tsv(repl_file, show_col_types = FALSE) %>%
    filter(tested) %>%
    summarise(pqtl_replication = paste(sprintf("%s, %s: %s", outcome, study,
                                               if_else(replicated, "replicated", "not replicated")),
                                       collapse = "; "), .by = protein) %>%
    rename(gene = protein)
} else {
  message("no pQTL replication (run 05_mr/05_pqtl_replication.R first)")
  tibble(gene = character(), pqtl_replication = character())
}

# Expression ----
# GTEx v10 median TPM: highest artery and highest tissue
gct <- fread(cmd = paste("gzip -dc", shQuote(p$gtex_median_tpm)), skip = 2, header = TRUE)
setnames(gct, 1:2, c("gene_id", "gtex_symbol"))
gct[, ensembl_id := sub("\\..*", "", gene_id)]
gct <- gct[!duplicated(ensembl_id)]
tissues <- setdiff(names(gct), c("gene_id", "gtex_symbol", "ensembl_id"))
arteries <- grep("^Artery", tissues, value = TRUE)
if (length(arteries) != 3) stop("expected three artery columns in the GTEx file, found: ", paste(arteries, collapse = ", "))
m <- as.matrix(gct[, ..tissues])
tissue_label <- function(x) str_squish(str_replace_all(sub("^Artery[ _-]+", "", x), "[_-]+", " "))
gtex <- tibble(ensembl_id = gct$ensembl_id, gtex_symbol = gct$gtex_symbol,
               artery_tpm = apply(m[, arteries], 1, max),
               artery_max = tissue_label(arteries[max.col(m[, arteries], ties.method = "first")]),
               top_tissue = str_squish(str_replace_all(tissues[max.col(m, ties.method = "first")], "[_]+", " ")),
               top_tpm = apply(m, 1, max))

# HCA arterial atlas cell type of highest specificity
hca_file <- file.path(cfg$results_dir, "vascular/hca_celltype_profile.tsv")
hca <- if (file.exists(hca_file)) {
  read_tsv(hca_file, show_col_types = FALSE) %>%
    filter(any(pct >= cfg$vascular$expressed_pct), .by = gene) %>%
    slice_max(specificity, n = 1, with_ties = FALSE, by = gene) %>%
    transmute(gene, hca_cell_type = cell_type, hca_specificity = specificity)
} else {
  message("no HCA profile (run 11_vascular_expression/03 first); HCA columns left empty")
  tibble(gene = character(), hca_cell_type = character(), hca_specificity = numeric())
}

# Table ----
tiers <- ev %>%
  select(gene, n_arteriopathies = n_traits, arteriopathies, in_table1) %>%
  left_join(gene_mr, by = "gene") %>%
  mutate(tier = coalesce(tier, 3L)) %>%
  left_join(disagree, by = "gene") %>%
  left_join(mono %>% select(gene, monogenic_disease, lof_mechanism), by = "gene") %>%
  mutate(mendelian_direction = case_when(
    !coalesce(lof_mechanism, FALSE) | is.na(implied_action) ~ NA_character_,
    implied_action == "Activate" ~ "Concordant",
    implied_action == "Inhibit" ~ "Discordant",
    TRUE ~ "Opposing across arteriopathies")) %>%
  left_join(tract %>% select(gene, ensembl_id, tractable, loeuf, burden_traits, burden_min_p, n_safety_events),
            by = "gene") %>%
  left_join(drug_status, by = "gene") %>%
  left_join(distinct(gtypes, ensembl_id, .keep_all = TRUE) %>% select(ensembl_id, gene_type), by = "ensembl_id") %>%
  mutate(gene_type = coalesce(gene_type, gtypes$gene_type[match(gene, gtypes$gene)])) %>%
  left_join(distinct(gtypes, ensembl_id, .keep_all = TRUE) %>% select(ensembl_id, gene_chr = chr, gene_start = start),
            by = "ensembl_id") %>%
  mutate(feasible_modalities = map2_chr(tractable, implied_action, function(t, a) {
    if (is.na(t) || is.na(a) || a == "Opposing") return(NA_character_)
    na_if(paste(intersect(str_split_1(t, ", "), p$modalities[[a]]), collapse = ", "), "")
  }),
  lof_intolerant = loeuf < p$loeuf_intolerant) %>%
  left_join(drug_match, by = "gene") %>%
  left_join(pav, by = "gene") %>%
  left_join(repl, by = "gene") %>%
  left_join(gtex %>% select(-gtex_symbol), by = "ensembl_id") %>%
  left_join(hca, by = "gene") %>%
  select(gene, tier, n_arteriopathies, arteriopathies, in_table1, mr_traits, coloc_traits,
         risk_up, risk_down, implied_action, direction_note, or_text, qtl, mr_effects, exposures_disagree,
         pqtl_protein_altering, pqtl_replication, gene_type, gene_chr, gene_start,
         monogenic_disease, mendelian_direction, tractable, feasible_modalities,
         max_clinical_phase, approved_drugs, matching_drugs, matching_max_phase,
         opposite_drugs, opposite_max_phase, loeuf, lof_intolerant, burden_traits, burden_min_p,
         n_safety_events, artery_tpm, artery_max, top_tissue, top_tpm, hca_cell_type,
         hca_specificity, ensembl_id) %>%
  arrange(tier, desc(n_arteriopathies), desc(coalesce(matching_max_phase, -1)), gene)

# Tier 1 loci: genes on the same chromosome with starts within locus_kb (single linkage)
locus_ids <- function(chr, pos, gap) {
  o <- order(chr, pos)
  id <- integer(length(pos))
  k <- 0L
  for (j in seq_along(o)) {
    i <- o[j]
    if (j == 1 || chr[i] != chr[o[j - 1]] || pos[i] - pos[o[j - 1]] > gap) k <- k + 1L
    id[i] <- k
  }
  id
}
t1_loci <- tiers %>%
  filter(tier == 1, !is.na(gene_start)) %>%
  mutate(locus_id = locus_ids(gene_chr, gene_start, p$locus_kb * 1e3)) %>%
  mutate(tier1_locus = sprintf("chr%s:%.1f Mb", first(gene_chr), min(gene_start) / 1e6),
         tier1_genes_at_locus = paste(sort(gene), collapse = ", "),
         n_tier1_at_locus = n(), .by = locus_id) %>%
  select(gene, tier1_locus, tier1_genes_at_locus, n_tier1_at_locus)
tiers <- tiers %>% left_join(t1_loci, by = "gene") %>% relocate(tier1_locus, n_tier1_at_locus, tier1_genes_at_locus, .after = tier)
loci_summary <- distinct(t1_loci, tier1_locus, n_tier1_at_locus, tier1_genes_at_locus)
message(nrow(loci_summary), " Tier 1 loci; ", sum(loci_summary$n_tier1_at_locus > 1), " with more than one Tier 1 gene")
print(loci_summary %>% filter(n_tier1_at_locus > 1) %>% arrange(desc(n_tier1_at_locus)), n = Inf)

missing_gtex <- sum(is.na(tiers$artery_tpm))
if (missing_gtex) message(missing_gtex, " genes without a GTEx median TPM (no Ensembl ID match)")
write_tsv(tiers, file.path(out, "target_tiers.tsv"))

# Main-text table: Tier 1 protein-coding genes with drugs acting on the target
phase_label <- function(x) case_when(x == 4 ~ "approved", x >= 1 ~ paste("phase", x))
salts <- c("ANHYDROUS", "SODIUM", "HYCLATE", "ESYLATE", "PHOSPHATE", "MESYLATE", "HYDROCHLORIDE",
           "DIHYDROCHLORIDE", "CALCIUM", "POTASSIUM", "ACETATE", "MALEATE", "TARTRATE", "SULFATE",
           "CITRATE", "SUCCINATE", "BESYLATE", "TOSYLATE", "FUMARATE", "HYDRATE", "MONOHYDRATE")
clean_drug <- function(x) {
  x <- str_squish(str_remove_all(toupper(x), paste0("\\b(", paste(salts, collapse = "|"), ")\\b")))
  if_else(str_detect(x, "\\d"), x, tolower(x))
}
best_drugs <- function(name, phase, n = 3) {
  if (length(name) == 0 || all(is.na(phase))) return(NA_character_)
  top <- max(phase, na.rm = TRUE)
  u <- sort(unique(clean_drug(name[which(phase == top)])))
  paste0(paste(head(u, n), collapse = ", "), if (length(u) > n) paste0(" +", length(u) - n) else "",
         " (", phase_label(top), ")")
}
main_drugs <- drugs %>%
  inner_join(tiers %>% filter(implied_action %in% c("Inhibit", "Activate")) %>% select(gene, implied_action),
             by = "gene") %>%
  filter(!is.na(drug_action)) %>%
  mutate(matches = drug_action == implied_action) %>%
  summarise(same = best_drugs(drug_name[matches], clinical_phase[matches]),
            opposite = best_drugs(drug_name[!matches], clinical_phase[!matches]),
            .by = gene)

main <- tiers %>%
  filter(tier == 1, gene_type == "protein_coding") %>%
  inner_join(main_drugs, by = "gene") %>%
  transmute(Gene = gene,
            `Arteriopathy (OR)` = or_text,
            QTL = qtl,
            `Implied action` = implied_action,
            `Drugs, same action` = same,
            `Drugs, opposite action` = opposite,
            `Artery TPM` = round(artery_tpm, 1),
            `Highest tissue (TPM)` = if_else(is.na(top_tissue), NA_character_,
                                             sprintf("%s (%s)", top_tissue, format(round(top_tpm), big.mark = ",", trim = TRUE))))
write_tsv(main, file.path(out, "target_tiers_main.tsv"))
notes <- tiers %>% filter(gene %in% main$Gene, !is.na(direction_note))
if (nrow(notes)) message("Table notes:\n", paste0(notes$gene, ": ", notes$direction_note, collapse = "\n"))

print(count(tiers, tier, implied_action))
print(tiers %>% filter(!is.na(mendelian_direction)) %>%
        select(gene, tier, monogenic_disease, mr_effects, mendelian_direction), n = Inf)
print(count(tiers, tier, gene_type))
print(main, n = Inf, width = Inf)

t1 <- filter(tiers, tier == 1)
write_counts("prioritization", c(
  genes = nrow(tiers),
  tier1 = nrow(t1), tier2 = sum(tiers$tier == 2), tier3 = sum(tiers$tier == 3),
  tier1_protein_coding = sum(t1$gene_type == "protein_coding", na.rm = TRUE),
  tier1_loci = nrow(loci_summary),
  tier1_loci_multigene = sum(loci_summary$n_tier1_at_locus > 1),
  tier1_pqtl_protein_altering = sum(!is.na(t1$pqtl_protein_altering)),
  tier1_2plus_arteriopathies = sum(t1$n_arteriopathies >= 2),
  main_table_genes = nrow(main),
  tier1_inhibit = sum(t1$implied_action == "Inhibit"),
  tier1_activate = sum(t1$implied_action == "Activate"),
  tier1_opposing = sum(t1$implied_action == "Opposing"),
  tier1_feasible = sum(!is.na(t1$feasible_modalities)),
  tier1_matching_drug = sum(!is.na(t1$matching_drugs)),
  tier1_matching_approved = sum(coalesce(t1$matching_max_phase == 4, FALSE)),
  tier1_opposite_drug = sum(!is.na(t1$opposite_drugs)),
  tier1_artery_tpm1 = sum(coalesce(t1$artery_tpm >= 1, FALSE)),
  tier1_lof_intolerant = sum(coalesce(t1$lof_intolerant, FALSE)),
  tier1_burden = sum(!is.na(t1$burden_traits))
))
