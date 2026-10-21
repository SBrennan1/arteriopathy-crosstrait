# Transcriptome-wide MR of GTEx v10 cis-eQTL instruments (aorta, coronary and tibial
# artery, whole blood, fibroblasts) on the five arteriopathies.

source("R/setup.R")
out <- res_dir("mr")

read_eqtl <- function(path) {
  d <- vroom(path, show_col_types = FALSE,
             col_types = cols(rsid = "c", effect_allele = "c", other_allele = "c"))
  gene_col <- intersect(c("exposure", "gene", "gene_symbol", "symbol", "Gene", "gene_name"), names(d))[1]
  if (is.na(gene_col)) stop("no gene column in ", path)
  if (!"gene_id" %in% names(d)) d$gene_id <- NA_character_
  d %>%
    transmute(exposure = .data[[gene_col]], gene_id = as.character(gene_id),
              rsid = as.character(rsid), beta, se, effect_allele, other_allele, pval,
              across(any_of("eaf")))
}

ivs <- map(set_names(names(cfg$mr$eqtl$qcd)),
           ~ read_eqtl(file.path(cfg$mr$eqtl$instruments, paste0(.x, ".tsv"))))

res <- map_dfr(traits, function(t) {
  outcome <- read_mr_outcome(t, rsids = unique(unlist(map(ivs, "rsid"))))
  map_dfr(ivs, run_mr, outcome = outcome, outcome_name = t, .id = "tissue")
})

gene_ids <- map_dfr(ivs, ~ .x %>% filter(!is.na(gene_id)) %>% distinct(exposure, .keep_all = TRUE) %>%
                      select(exposure, gene_id), .id = "tissue")

res %>%
  left_join(gene_ids, by = c("tissue", "exposure")) %>%
  relocate(gene_id, .after = exposure) %>%
  write_tsv(file.path(out, "mr_eqtl_unadjusted.tsv"))
