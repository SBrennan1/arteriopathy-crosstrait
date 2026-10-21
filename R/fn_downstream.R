# Helpers for the downstream analyses

read_stage <- function(...) {
  read_tsv(file.path(cfg$results_dir, ...), show_col_types = FALSE, guess_max = 1e5)
}

# Gene-outcome pairs with FDR-significant MR: pQTL on either platform, eQTL in arterial
# tissues
mr_hits <- function(fdr = cfg$mr$fdr) {
  list(
    pqtl = read_stage("mr", "mr_pqtl.tsv") %>%
      filter(qval < fdr) %>%
      distinct(gene = exposure, trait = outcome),
    eqtl = read_stage("mr", "mr_eqtl.tsv") %>%
      filter(tissue %in% cfg$mr$eqtl$arterial_tissues, qval < fdr) %>%
      distinct(gene = exposure, trait = outcome)
  )
}

# Gene-level evidence across ASSET loci, pQTL MR and arterial eQTL MR, and the
# arteriopathies associated with each gene
convergence_evidence <- function(pqtl, eqtl, p = cfg$convergence) {
  asset <- read_stage("v2g", "asset_v2g.tsv") %>% filter(!is.na(gene))
  pqtl_coloc <- read_stage("qtl_coloc", "pqtl_coloc.tsv")
  eqtl_coloc <- read_stage("qtl_coloc", "eqtl_coloc.tsv")

  # Optional exclusions (config: asset_trait_exclude)
  excl <- p$asset_trait_exclude %||% list()
  exclude <- tibble(gene = as.character(rep(names(excl), lengths(excl))),
                    trait = as.character(unlist(excl, use.names = FALSE)))
  asset_traits <- asset %>%
    transmute(gene, trait = traits) %>%
    separate_rows(trait, sep = ",\\s*") %>%
    filter(trait != "") %>%
    distinct() %>%
    anti_join(exclude, by = c("gene", "trait"))

  traits_by_gene <- bind_rows(asset_traits, pqtl, eqtl) %>%
    distinct(gene, trait) %>%
    group_by(gene) %>%
    summarise(n_traits = n_distinct(trait),
              arteriopathies = paste(sort(unique(trait)), collapse = ", "), .groups = "drop")

  coloc_p <- pqtl_coloc %>%
    filter(coalesce(abf_pph4 >= p$coloc_pph4, FALSE) | coalesce(susie_pph4 >= p$coloc_pph4, FALSE)) %>%
    pull(protein)
  coloc_e <- eqtl_coloc %>% filter(PP.H4 >= p$coloc_pph4) %>% pull(gene)

  genes <- tibble(gene = unique(c(asset$gene, pqtl$gene, eqtl$gene))) %>%
    mutate(ASSET = gene %in% asset$gene,
           pQTL_MR = gene %in% pqtl$gene,
           pQTL_coloc = gene %in% coloc_p,
           eQTL_MR = gene %in% eqtl$gene,
           eQTL_coloc = gene %in% coloc_e,
           n_modalities = ASSET + pQTL_MR + eQTL_MR) %>%
    left_join(traits_by_gene, by = "gene") %>%
    mutate(in_table1 = coalesce(n_traits >= p$min_traits, FALSE) & !str_detect(gene, "^ENSG"))

  list(genes = genes, traits = bind_rows(asset_traits %>% mutate(source = "ASSET"),
                                         pqtl %>% mutate(source = "pQTL_MR"),
                                         eqtl %>% mutate(source = "eQTL_MR")))
}

# Biological theme of a pathway / gene-set name (first matching rule wins)
pathway_theme <- function(description) {
  rules <- c(
    "ECM/Structural" = "collagen|laminin|elastin|extracellular matrix|extracellular structure|external encapsulating|ecm|integrin|fibril|elastic fibre|fibronectin|ossification",
    "TGF-β/BMP Signaling" = "tgf|smad|bmp|activin",
    "Complement System" = "complement|c4 and c2|ficolin|lectin pathway",
    "Growth Factor Signaling" = "pdgf|vegf|fgf|egf|growth factor|igf",
    "Lipid/Scavenger" = "scaveng|lipoprotein|lipid|cholesterol|sterol|triglyceride|acylglycerol|phospholipid|fatty acid|ppar",
    "Inflammation/Immune" = "inflam|cytokine|interleukin|immune|interferon|mhc|antigen process|efferocytosis",
    "Cell Cycle/DNA Damage" = "cell cycle|mitotic|checkpoint|cyclin|cdk|dna damage|tp53|chk1|chk2|telomere",
    "Vascular" = "vascular|vasculature|angiogenesis|blood vessel|endothel|artery morpho|circulatory system|blood.brain barrier|nitric oxide|guanylate cyclase",
    "Smooth Muscle" = "smooth muscle",
    "Cell Migration" = "migration|motility|locomotion|chemotaxis",
    "Cell Adhesion" = "adhesion|focal adhesion|gap junction",
    "Cytoskeleton" = "actin|myosin|contraction|cytoskelet",
    "Semaphorin/Guidance" = "semaphorin|axon guid",
    "Endocytosis/Vesicular Transport" = "endocytosis|vesicle.mediated transport|receptor.mediated endocytosis|import into cell|regulation of endocytosis",
    "Platelet/Hemostasis" = "platelet|coagulation|hemostasis|thrombosis",
    "Glycosylation" = "glycosylat|b3galtl|tsr domain",
    "Ion Transport" = "sodium ion transport|cation transport|slc.mediated|ion channel",
    "RAAS/Blood Pressure" = "renin|aldosterone|angiotensin",
    "Wnt Signaling" = "wnt|ctnnb|beta-catenin|tcf|lef",
    "Hedgehog Signaling" = "hedgehog|gli |hh responsive",
    "Cell Proliferation" = "proliferation",
    "Signal Transduction" = "signal transduction|signaling|signalling|cell communication",
    "Proteolysis/Remodeling" = "proteolysis|protease|peptidase|degradation|remodeling",
    "Amyloid Processing" = "amyloid",
    "Development" = "kidney develop|renal system|embryo|morphogenesis|endocrine system develop"
  )
  first_match(tolower(description), rules, names(rules), default = "Other")
}

# Label of the first pattern matching each element of x
first_match <- function(x, patterns, labels, default = NA_character_, ignore_case = FALSE) {
  out <- rep(NA_character_, length(x))
  for (i in seq_along(patterns)) {
    hit <- is.na(out) & grepl(patterns[i], x, ignore.case = ignore_case)
    out[hit] <- labels[i]
  }
  out[is.na(out)] <- default
  out
}

# Cytogenetic band at a GRCh37 position (Ensembl REST overlap endpoint)
ensembl_cytoband <- function(chr, pos, url = cfg$annotation$ensembl_rest) {
  r <- httr::GET(sprintf("%s/overlap/region/human/%s:%d-%d", url, chr, as.integer(pos), as.integer(pos)),
                 query = list(feature = "band", `content-type` = "application/json"))
  Sys.sleep(0.2)
  if (httr::status_code(r) != 200) return(NA_character_)
  band <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))
  if (length(band) == 0 || is.null(band$id)) NA_character_ else paste0(chr, band$id[1])
}
