# coloc.abf between GTEx v10 eQTLs and the outcome GWAS for each FDR-significant eQTL MR
# association (gene, tissue and outcome).

source("R/setup.R")
p <- cfg$qtl_coloc
out <- res_dir("qtl_coloc")

hits <- read_tsv(file.path(cfg$results_dir, "mr/mr_eqtl.tsv"), show_col_types = FALSE) %>%
  filter(qval < p$fdr) %>%
  distinct(outcome, tissue, gene = exposure, gene_id)

eqtl_files <- list.files(p$eqtl$dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
chain <- rtracklayer::import.chain(cfg$chain_hg38_to_hg19)

res <- map_dfr(unique(hits$outcome), function(t) {
  gwas <- read_univar(t, c("chr", "pos", "effect_allele", "other_allele", "beta", "se")) %>%
    as.data.table()
  ncase <- cfg$traits[[t]]$ncase
  ncontrol <- cfg$traits[[t]]$ncontrol

  pmap_dfr(filter(hits, outcome == t), function(outcome, tissue, gene, gene_id) {
    r <- tryCatch(
      read_eqtl_subset(tissue, gene_id, eqtl_files, chain) %>%
        coloc_eqtl(gwas, ncase, ncontrol) %>%
        mutate(reason = NA_character_),
      error = function(e) tibble(reason = conditionMessage(e)))
    bind_cols(tibble(outcome, tissue, gene, gene_id), r)
  })
})

write_tsv(res %>% filter(!is.na(reason)) %>% select(outcome, tissue, gene, gene_id, reason),
          file.path(out, "eqtl_coloc_not_run.tsv"))

write_tsv(res %>%
            filter(is.na(reason)) %>%
            select(outcome, tissue, gene, gene_id, nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4,
                   method) %>%
            arrange(outcome, tissue, gene),
          file.path(out, "eqtl_coloc.tsv"))
