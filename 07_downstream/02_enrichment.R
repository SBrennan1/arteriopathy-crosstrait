# Pathway enrichment (GO BP, KEGG, Reactome; clusterProfiler) of the FDR-significant
# pQTL MR proteins and arterial eQTL MR genes, pooled across outcomes, against all
# proteins or arterial genes tested.

source("R/setup.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})
select <- dplyr::select

out <- res_dir("enrichment")
p <- cfg$enrichment
fdr <- cfg$mr$fdr

pqtl <- read_stage("mr", "mr_pqtl.tsv")
eqtl <- read_stage("mr", "mr_eqtl.tsv") %>% filter(tissue %in% cfg$mr$eqtl$arterial_tissues)

to_entrez <- function(symbols) {
  if (length(symbols) == 0) return(character())
  ids <- suppressWarnings(suppressMessages(
    bitr(unique(symbols), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)))
  unique(ids$ENTREZID)
}

universe <- list(pQTL_MR = to_entrez(pqtl$exposure), Arterial_eQTL_MR = to_entrez(eqtl$exposure))

enrich <- function(symbols, source) {
  genes <- to_entrez(symbols)
  bg <- universe[[source]]
  args <- list(gene = genes, universe = bg, pAdjustMethod = "BH", pvalueCutoff = p$padj,
               qvalueCutoff = p$qvalue_cutoff)
  fits <- compact(list(
    GO_BP = tryCatch(do.call(enrichGO, c(args, list(OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE))),
                     error = function(e) NULL),
    KEGG = tryCatch(do.call(enrichKEGG, c(args, list(organism = "hsa"))), error = function(e) NULL),
    Reactome = tryCatch(do.call(enrichPathway, c(args, list(organism = "human", readable = TRUE))),
                        error = function(e) NULL)
  ))
  imap_dfr(fits, ~ as_tibble(.x@result) %>% mutate(database = .y)) %>%
    filter(p.adjust < p$padj) %>%
    transmute(source = !!source, database, id = ID, description = Description,
              gene_ratio = GeneRatio, bg_ratio = BgRatio, pvalue, p_adjust = p.adjust, qvalue,
              genes = geneID, count = Count, theme = pathway_theme(description))
}

gene_lists <- list(pQTL_MR = unique(pqtl$exposure[pqtl$qval < fdr]),
                   Arterial_eQTL_MR = unique(eqtl$exposure[eqtl$qval < fdr]))

imap_dfr(gene_lists, ~ tibble(source = .y, n_genes = length(.x),
                              genes = paste(sort(.x), collapse = ", "))) %>%
  write_tsv(file.path(out, "mr_gene_lists.tsv"))

terms <- imap_dfr(gene_lists, enrich) %>% arrange(source, database, p_adjust)
write_tsv(terms, file.path(out, "mr_pathway_enrichment.tsv"))

write_counts("enrichment_mr", c(
  pqtl_genes = length(gene_lists$pQTL_MR),
  arterial_eqtl_genes = length(gene_lists$Arterial_eQTL_MR),
  pqtl_terms = sum(terms$source == "pQTL_MR"),
  arterial_eqtl_terms = sum(terms$source == "Arterial_eQTL_MR")
))
