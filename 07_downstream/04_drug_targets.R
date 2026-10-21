# Open Targets Platform annotation of genes implicated by ASSET, MTAG, pQTL MR or eQTL MR:
# tractability, drugs and clinical candidates, gnomAD constraint, safety and rare-variant
# burden evidence. Responses are cached (ot_raw.rds); set OT_REFRESH=1 to requery.

source("R/setup.R")
out <- res_dir("drug_targets")
ot_url <- cfg$v2g$opentargets_url
p <- cfg$drug_targets

hits <- mr_hits()
ev <- convergence_evidence(hits$pqtl, hits$eqtl)

pipeline <- ev$genes %>%
  select(gene, ASSET, pQTL_MR, eQTL_MR) %>%
  left_join(ev$traits %>%
              summarise(pipeline_traits = paste(sort(unique(trait)), collapse = "; "), .by = gene),
            by = "gene")

mtag <- read_stage("mtag", "mtag_leads.tsv") %>%
  inner_join(read_stage("v2g", "mtag_v2g.tsv") %>% select(rsid, gene, evidence), by = "rsid") %>%
  filter(!is.na(gene)) %>%
  summarise(mtag_traits = paste(sort(unique(target_trait)), collapse = "; "),
            mtag_pairs = paste(sort(unique(mtag_pair)), collapse = "; "),
            mtag_mapping = paste(sort(unique(evidence)), collapse = "; "),
            mtag_novel = any(novel, na.rm = TRUE), .by = gene)

join_traits <- function(a, b) {
  t <- unique(unlist(str_split(c(a, b), ";\\s*")))
  paste(sort(t[t != ""]), collapse = "; ")
}

master <- tibble(gene = unique(c(pipeline$gene, mtag$gene))) %>%
  left_join(pipeline, by = "gene") %>%
  left_join(mtag, by = "gene") %>%
  mutate(in_asset = coalesce(ASSET, FALSE),
         in_pqtl = coalesce(pQTL_MR, FALSE),
         in_eqtl = coalesce(eQTL_MR, FALSE),
         in_mtag = gene %in% mtag$gene,
         associated_traits = map2_chr(coalesce(pipeline_traits, ""), coalesce(mtag_traits, ""), join_traits),
         n_traits = map_int(str_split(associated_traits, ";\\s*"), ~ sum(.x != "")),
         n_modalities = in_asset + in_pqtl + in_eqtl + in_mtag,
         discovery_sources = pmap_chr(list(in_asset, in_mtag, in_pqtl, in_eqtl), function(a, m, p, e) {
           paste(c("ASSET", "MTAG", "pQTL", "eQTL")[c(a, m, p, e)], collapse = "; ")
         }),
         mtag_only = in_mtag & !in_asset & !in_pqtl & !in_eqtl,
         v2g_evidence = mtag_mapping) %>%
  filter(!str_detect(gene, "^ENSG")) %>%
  arrange(desc(n_traits), desc(n_modalities), gene)

write_tsv(master %>%
            select(gene, n_traits, n_modalities, discovery_sources, in_asset, in_mtag, in_pqtl,
                   in_eqtl, mtag_only, mtag_novel, associated_traits, v2g_evidence),
          file.path(out, "master_gene_list.tsv"))

# Open Targets queries ----
ot_query <- function(query, tries = 3) {
  for (i in seq_len(tries)) {
    r <- tryCatch(httr::POST(ot_url, body = list(query = query), encode = "json", httr::timeout(60)),
                  error = function(e) NULL)
    if (!is.null(r)) {
      x <- tryCatch(httr::content(r, as = "parsed", simplifyVector = FALSE), error = function(e) NULL)
      if (length(x$errors)) stop("Open Targets: ", x$errors[[1]]$message, call. = FALSE)
      if (httr::status_code(r) == 200) return(x$data)
    }
    Sys.sleep(2 * i)
  }
  stop("Open Targets API not reachable", call. = FALSE)
}
s_ <- function(x) if (length(x) == 0) NA_character_ else paste(unlist(x), collapse = "; ")
n_ <- function(x) if (length(x) == 0) NA_real_ else as.numeric(x[[1]])
collapse_ <- function(x) na_if(paste(unique(na.omit(x[x != ""])), collapse = "; "), "")

# Ensembl IDs from GENCODE by symbol, otherwise the Open Targets search
gencode <- fread(cmd = paste("gzip -dc", shQuote(cfg$gtex_v10$gencode), "| awk -F'\\t' '$3==\"gene\"'"),
                 header = FALSE, sep = "\t", select = 9, col.names = "attr")[!grepl("_PAR_Y", attr)]
sym <- gencode[, .(gene = sub('.*gene_name "([^"]+)".*', "\\1", attr),
                   ensembl_id = sub('.*gene_id "(ENSG\\d+).*', "\\1", attr))]
sym <- sym[!gene %in% gene[duplicated(gene)]]
ot_search <- function(symbol) {
  h <- ot_query(sprintf('{ search(queryString: "%s", entityNames: ["target"], page: {size: 1, index: 0}) {
    hits { id } } }', symbol))$search$hits
  if (length(h)) h[[1]]$id else NA_character_
}
ids <- tibble(gene = master$gene) %>%
  left_join(as_tibble(sym), by = "gene") %>%
  mutate(id_source = if_else(is.na(ensembl_id), "Open Targets search", "GENCODE v39"))
miss <- is.na(ids$ensembl_id)
ids$ensembl_id[miss] <- map_chr(ids$gene[miss], ~ tryCatch(ot_search(.x), error = function(e) NA_character_))
message(sum(!miss), " genes matched in GENCODE v39, ", sum(miss), " searched, ",
        sum(is.na(ids$ensembl_id)), " without an Ensembl ID")

# Diseases for the burden look-up, with their descendants
dis <- imap_dfr(p$diseases, function(d, trait) {
  id <- d$id %||% ot_query(sprintf('{ search(queryString: "%s", entityNames: ["disease"], page: {size: 1, index: 0}) {
    hits { id } } }', d$query))$search$hits[[1]]$id
  x <- ot_query(sprintf('{ disease(efoId: "%s") { id name descendants } }', id))$disease
  tibble(trait, query = d$query %||% NA_character_, efo_id = x$id, efo_name = x$name,
         n_descendants = length(x$descendants), ids = list(c(x$id, unlist(x$descendants))))
})
print(select(dis, -ids))
write_tsv(select(dis, -ids), file.path(out, "ot_diseases.tsv"))
efo <- paste0('"', unique(unlist(dis$ids)), '"', collapse = ", ")

target_query <- '{ target(ensemblId: "%s") {
  id approvedSymbol isEssential
  tractability { label modality value }
  prioritisation { items { key value } }
  geneticConstraint { constraintType oeUpper upperBin6 }
  safetyLiabilities { event datasource }
  chemicalProbes { id isHighQuality }
  drugAndClinicalCandidates { count rows {
    maxClinicalStage
    diseases { diseaseFromSource disease { id name therapeuticAreas { id } } }
    drug { id name drugType maximumClinicalStage
      mechanismsOfAction { rows { mechanismOfAction actionType targets { id } } } } } }
  burden: evidences(efoIds: [%s], datasourceIds: ["gene_burden"], size: 100) { count rows {
    disease { id name } projectId cohortId studyId statisticalMethod ancestry beta oddsRatio
    pValueMantissa pValueExponent directionOnTrait studyCasesWithQualifyingVariants } }
  associatedDiseases(page: {index: 0, size: 10}) { rows { score disease { name } } }
} }'

cache_file <- file.path(out, "ot_raw.rds")
raw <- if (file.exists(cache_file) && Sys.getenv("OT_REFRESH") == "") readRDS(cache_file) else list()
todo <- setdiff(unique(na.omit(ids$ensembl_id)), names(raw))
if (length(todo) > 0) {
  # Genes not in Open Targets are stored as FALSE; failed queries are retried on rerun
  get_target <- function(id) ot_query(sprintf(target_query, id, efo))$target %||% FALSE
  raw[[todo[1]]] <- get_target(todo[1])
  for (id in todo[-1]) {
    raw[id] <- list(tryCatch(get_target(id), error = function(e) { message(id, ": ", conditionMessage(e)); NULL }))
    Sys.sleep(0.2)
  }
  raw <- compact(raw)
  version <- ot_query("{ meta { apiVersion { x y z } dataVersion { year month } } }")$meta
  attr(raw, "version") <- tibble(api = paste(unlist(version$apiVersion), collapse = "."),
                                 data = paste(unlist(version$dataVersion), collapse = "."),
                                 queried = as.character(Sys.Date()))
  saveRDS(raw, cache_file)
}
write_tsv(attr(raw, "version"), file.path(out, "opentargets_version.tsv"))
failed <- setdiff(unique(na.omit(ids$ensembl_id)), names(raw))
if (length(failed)) message(length(failed), " queries failed; rerun to retry them")
not_in_ot <- ids %>% filter(ensembl_id %in% names(raw)[map_lgl(raw, isFALSE)])
if (nrow(not_in_ot)) message("not in Open Targets: ", paste(not_in_ot$gene, collapse = ", "))
raw <- keep(raw, is.list)

# Parse ----
tract <- map_dfr(raw, function(t) map_dfr(t$tractability, ~ tibble(
  ensembl_id = t$id, label = s_(.x$label), modality = s_(.x$modality), value = isTRUE(.x$value))))

# Clinical stage to phase (approved = 4; combined phases take the higher)
stage_phase <- function(s) {
  num <- map_dbl(str_extract_all(coalesce(s, ""), "\\d"), ~ if (length(.x)) max(as.numeric(.x)) else NA_real_)
  if_else(str_detect(coalesce(s, ""), "APPROV"), 4, num)
}
# Cardiovascular and lipid indication flags
cv_area <- p$cardiovascular_area
flag <- function(x, cat) coalesce(str_detect(tolower(x), paste(p$indication_patterns[[cat]], collapse = "|")), FALSE)
drugs <- map_dfr(raw, function(t) map_dfr(t$drugAndClinicalCandidates$rows, function(r) {
  moa <- keep(r$drug$mechanismsOfAction$rows %||% list(),
              ~ t$id %in% map_chr(.x$targets %||% list(), ~ s_(.x$id)))
  dz <- r$diseases %||% list()
  areas <- unlist(map(dz, ~ map_chr(.x$disease$therapeuticAreas %||% list(), ~ s_(.x$id))))
  tibble(ensembl_id = t$id, drug_id = s_(r$drug$id), drug_name = s_(r$drug$name),
         drug_type = s_(r$drug$drugType), clinical_stage = s_(r$drug$maximumClinicalStage),
         mechanism = collapse_(map_chr(moa, ~ s_(.x$mechanismOfAction))),
         action_type = collapse_(map_chr(moa, ~ s_(.x$actionType))),
         indications = collapse_(map_chr(dz, ~ s_(.x$disease$name %||% .x$diseaseFromSource))),
         cv_area = cv_area %in% areas)
})) %>%
  summarise(across(c(drug_name, drug_type, clinical_stage), first),
            across(c(mechanism, action_type, indications), ~ collapse_(unlist(str_split(.x, "; ")))),
            cv_area = any(cv_area), .by = c(ensembl_id, drug_id)) %>%
  mutate(max_phase = stage_phase(clinical_stage),
         cv_indication = cv_area | flag(indications, "cardiovascular"),
         lipid_indication = flag(indications, "lipid"))
if (nrow(drugs) == 0) stop("no drugs returned for any gene; check the drugAndClinicalCandidates query")
print(count(drugs, clinical_stage, max_phase))

prio <- map_dfr(raw, function(t) map_dfr(t$prioritisation$items %||% list(), ~ tibble(
  ensembl_id = t$id, key = s_(.x$key), value = s_(.x$value))))
if (nrow(prio)) prio <- pivot_wider(prio, names_from = key, values_from = value, values_fn = first)

info <- map_dfr(raw, function(t) {
  lof <- keep(t$geneticConstraint %||% list(), ~ identical(.x$constraintType, "lof"))
  tibble(ensembl_id = t$id, ot_symbol = s_(t$approvedSymbol),
         depmap_essential = isTRUE(t$isEssential),
         loeuf = if (length(lof)) n_(lof[[1]]$oeUpper) else NA_real_,
         loeuf_bin6 = if (length(lof)) n_(lof[[1]]$upperBin6) else NA_real_,
         n_safety_events = length(t$safetyLiabilities),
         safety_events = collapse_(map_chr(t$safetyLiabilities %||% list(), ~ s_(.x$event))),
         n_hq_chemical_probes = sum(map_lgl(t$chemicalProbes %||% list(), ~ isTRUE(.x$isHighQuality))),
         top_associated_diseases = collapse_(map_chr(t$associatedDiseases$rows %||% list(), ~ s_(.x$disease$name))))
})

trait_of <- function(id) collapse_(dis$trait[map_lgl(dis$ids, ~ id %in% .x)])
burden <- map_dfr(raw, function(t) map_dfr(t$burden$rows %||% list(), function(e) tibble(
  ensembl_id = t$id, disease_id = s_(e$disease$id), disease = s_(e$disease$name),
  project = s_(e$projectId), cohort = s_(e$cohortId), study_id = s_(e$studyId),
  method = s_(e$statisticalMethod), ancestry = s_(e$ancestry), beta = n_(e$beta),
  odds_ratio = n_(e$oddsRatio), p = n_(e$pValueMantissa) * 10^n_(e$pValueExponent),
  direction_on_trait = s_(e$directionOnTrait), cases_with_variants = n_(e$studyCasesWithQualifyingVariants))))
if (nrow(burden)) {
  burden <- mutate(burden, trait = map_chr(disease_id, trait_of))
} else {
  burden <- tibble(ensembl_id = character(), disease = character(), trait = character(), p = numeric())
}

# Genes with a different Open Targets symbol
renamed <- ids %>% inner_join(info, by = "ensembl_id") %>% filter(gene != ot_symbol)
if (nrow(renamed)) print(select(renamed, gene, ot_symbol, ensembl_id, id_source), n = Inf)

# Tables ----
genes <- ids %>% filter(!is.na(ensembl_id)) %>% select(gene, ensembl_id)

tract_wide <- tract %>%
  filter(value) %>%
  mutate(modality = recode(modality, SM = "tractability_sm", AB = "tractability_ab",
                           PR = "tractability_protac", OC = "tractability_oc",
                           .default = paste0("tractability_", tolower(modality)))) %>%
  distinct(ensembl_id, modality) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = modality, values_from = value, values_fill = FALSE)
for (col in c("tractability_sm", "tractability_ab", "tractability_protac", "tractability_oc")) {
  if (!col %in% names(tract_wide)) tract_wide[[col]] <- logical(nrow(tract_wide))
}
tract_labels <- tract %>%
  filter(value) %>%
  summarise(tractability_labels = paste(unique(label), collapse = "; "), .by = ensembl_id)

drug_summary <- drugs %>%
  summarise(n_drugs = n_distinct(drug_id),
            max_clinical_phase = suppressWarnings(max(max_phase, na.rm = TRUE)),
            approved_drugs = collapse_(drug_name[which(max_phase == 4)]),
            phase3_drugs = collapse_(drug_name[which(max_phase == 3)]),
            cv_drugs = collapse_(drug_name[cv_indication]),
            lipid_drugs = collapse_(drug_name[lipid_indication]),
            all_mechanisms = collapse_(mechanism), .by = ensembl_id) %>%
  mutate(max_clinical_phase = if_else(is.infinite(max_clinical_phase), NA_real_, max_clinical_phase))

burden_summary <- if (nrow(burden)) {
  burden %>%
    summarise(burden_traits = collapse_(unlist(str_split(trait, "; "))),
              burden_min_p = suppressWarnings(min(p, na.rm = TRUE)), .by = ensembl_id) %>%
    mutate(burden_min_p = if_else(is.infinite(burden_min_p), NA_real_, burden_min_p))
} else tibble(ensembl_id = character(), burden_traits = character(), burden_min_p = numeric())

master %>%
  select(gene, n_traits, associated_traits, discovery_sources, asset = in_asset, mtag = in_mtag,
         pqtl = in_pqtl, eqtl = in_eqtl, mtag_novel, v2g_evidence) %>%
  left_join(ids %>% select(gene, ensembl_id), by = "gene") %>%
  left_join(info, by = "ensembl_id") %>%
  left_join(tract_wide, by = "ensembl_id") %>%
  left_join(tract_labels, by = "ensembl_id") %>%
  left_join(drug_summary, by = "ensembl_id") %>%
  left_join(burden_summary, by = "ensembl_id") %>%
  mutate(across(starts_with("tractability_") & where(is.logical), ~ coalesce(.x, FALSE)),
         n_drugs = coalesce(n_drugs, 0L)) %>%
  select(gene, ensembl_id, ot_symbol, n_traits, associated_traits, discovery_sources, asset, mtag,
         pqtl, eqtl, mtag_novel, v2g_evidence, tractability_sm, tractability_ab, tractability_protac,
         tractability_oc, tractability_labels, n_drugs, max_clinical_phase, approved_drugs,
         phase3_drugs, cv_drugs, lipid_drugs, all_mechanisms, loeuf, loeuf_bin6, depmap_essential,
         n_safety_events, safety_events, n_hq_chemical_probes, burden_traits, burden_min_p,
         top_associated_diseases) %>%
  arrange(desc(n_traits), desc(n_drugs), gene) %>%
  write_tsv(file.path(out, "drug_tractability.tsv"))

drug_details <- genes %>%
  inner_join(drugs, by = "ensembl_id") %>%
  transmute(gene, drug_name, drug_id, drug_type, clinical_phase = max_phase, clinical_stage,
            mechanism_of_action = mechanism, action_type, indications, cv_indication, lipid_indication,
            status = if_else(clinical_phase == 4, "Approved", "Investigational", missing = NA_character_)) %>%
  arrange(gene, desc(clinical_phase), drug_name)
write_tsv(drug_details, file.path(out, "drug_details.tsv"))

if (nrow(prio)) write_tsv(genes %>% inner_join(prio, by = "ensembl_id"), file.path(out, "target_prioritisation.tsv"))
write_tsv(genes %>% inner_join(burden, by = "ensembl_id") %>% arrange(p), file.path(out, "burden_evidence.tsv"))

write_counts("drug_targets", c(
  genes_queried = nrow(master),
  genes_asset = sum(master$in_asset),
  genes_mtag = sum(master$in_mtag),
  genes_mtag_novel = sum(master$mtag_novel, na.rm = TRUE),
  genes_pqtl = sum(master$in_pqtl),
  genes_eqtl = sum(master$in_eqtl),
  genes_mtag_only = sum(master$mtag_only),
  genes_mapped_ensembl = sum(!is.na(ids$ensembl_id)),
  genes_failed = length(failed),
  genes_not_in_opentargets = nrow(not_in_ot),
  genes_with_drugs = n_distinct(drug_details$gene),
  genes_with_approved_drugs = n_distinct(drug_details$gene[which(drug_details$clinical_phase == 4)]),
  genes_with_cv_drugs = n_distinct(drug_details$gene[drug_details$cv_indication]),
  genes_with_lipid_drugs = n_distinct(drug_details$gene[drug_details$lipid_indication]),
  gene_drug_pairs = nrow(drug_details),
  genes_with_burden = n_distinct(burden$ensembl_id),
  genes_with_safety_events = sum(info$n_safety_events > 0)
))
