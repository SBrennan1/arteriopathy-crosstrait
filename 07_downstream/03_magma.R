# MAGMA gene and gene-set analysis of the ASSET cross-trait statistics, and FUMA
# GENE2FUNC tissue and cellular component enrichment of the ASSET genes (FUMA output).

source("R/setup.R")
out <- res_dir("enrichment")
ecfg <- cfg$enrichment

read_magma <- function(path) {
  d <- fread(path, skip = sum(grepl("^#", readLines(path))), header = TRUE, fill = TRUE)
  setnames(d, toupper(names(d)))
  if (!"FULL_NAME" %in% names(d) && "VARIABLE" %in% names(d)) d$FULL_NAME <- d$VARIABLE
  as_tibble(d)
}

# MAGMA ----
genes <- read_magma(file.path(ecfg$asset_fuma_dir, "magma.genes.out"))
gsa <- read_magma(file.path(ecfg$asset_fuma_dir, "magma.gsa.out"))

genes %>%
  mutate(p_bonferroni = pmin(P * n(), 1)) %>%
  filter(P < 0.05 / n()) %>%
  arrange(P) %>%
  transmute(symbol = SYMBOL, ensembl_id = GENE, chr = CHR, start = START, stop = STOP,
            nsnps = NSNPS, nparam = NPARAM, zstat = ZSTAT, p = P, p_bonferroni,
            p_adjust_method = sprintf("Bonferroni over %d genes", nrow(genes))) %>%
  write_tsv(file.path(out, "magma_asset_genes.tsv"))

msigdb_db <- function(x) {
  first_match(x, c("^GOBP_", "^GOMF_", "^GOCC_", "^KEGG_", "^REACTOME_", "^WP_"),
              c("GO Biological Process", "GO Molecular Function", "GO Cellular Component",
                "KEGG", "Reactome", "WikiPathways"), default = "Other")
}

# Canonical gene sets (GO, KEGG, Reactome, WikiPathways; 10-500 genes), top 200
canonical <- gsa %>%
  filter(str_detect(FULL_NAME, "^GO|^KEGG|^REACTOME|^WP_"),
         NGENES >= ecfg$asset_canonical_ngenes[1], NGENES <= ecfg$asset_canonical_ngenes[2]) %>%
  mutate(p_bonferroni = pmin(P * n(), 1), p_fdr = p.adjust(P, "BH")) %>%
  arrange(P) %>%
  transmute(pathway_id = FULL_NAME,
            pathway = str_to_title(str_replace_all(
              str_remove(FULL_NAME, "^(GOBP|GOMF|GOCC|KEGG|REACTOME|WP)_"), "_", " ")),
            database = msigdb_db(FULL_NAME), n_genes = NGENES, beta = BETA,
            beta_std = if ("BETA_STD" %in% names(gsa)) BETA_STD else NA_real_, se = SE, p = P,
            p_fdr, p_bonferroni)
write_tsv(head(canonical, 200), file.path(out, "magma_asset_pathways_canonical.tsv"))

# FUMA GENE2FUNC ----
g2f <- ecfg$gene2func_dir
gs <- read_tsv(file.path(g2f, "GS.txt"), show_col_types = FALSE)
deg <- read_tsv(file.path(g2f, "gtex_v8_ts_DEG.txt"), show_col_types = FALSE)

deg %>%
  filter(Category == "DEG.up") %>%
  arrange(adjP) %>%
  transmute(tissue = str_replace_all(GeneSet, "_", " "), n_genes_in_set = N_genes,
            n_overlapping = N_overlap, p, adj_p = adjP) %>%
  write_tsv(file.path(out, "gene2func_tissue_deg.tsv"))

gs %>%
  filter(Category == "GO_cc") %>%
  arrange(adjP) %>%
  transmute(gene_set = str_to_title(str_replace_all(str_remove(GeneSet, "^GOCC_"), "_", " ")),
            n_genes_in_set = N_genes, n_overlapping = N_overlap, p, adj_p = adjP,
            genes = str_replace_all(genes, ":", ", ")) %>%
  write_tsv(file.path(out, "gene2func_go_cc.tsv"))

write_counts("enrichment_magma", c(
  asset_magma_genes_bonferroni = sum(genes$P < 0.05 / nrow(genes)),
  asset_magma_canonical_bonferroni = sum(canonical$p_bonferroni < 0.05)
))
