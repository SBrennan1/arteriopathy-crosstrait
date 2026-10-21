# Flags FDR-significant pQTL instruments that are, or tag (r2 >= 0.8, 1000 Genomes EUR),
# a protein-altering variant in the measured protein (possible epitope effects).
# Ensembl GRCh37 REST (LD and VEP).

source("R/setup.R")
out <- res_dir("prioritization")
url <- cfg$annotation$ensembl_rest
p <- cfg$prioritization$protein_altering

inst <- read_stage("qtl_coloc", "pqtl_mr_coloc.tsv") %>%
  filter(qval < cfg$mr$fdr) %>%
  separate_rows(snps, sep = ";") %>%
  distinct(gene = exposure, rsid = snps) %>%
  filter(str_detect(rsid, "^rs"))
message(nrow(inst), " protein x instrument pairs, ", n_distinct(inst$rsid), " variants")

rest_get <- function(path, query = list()) {
  for (i in 1:3) {
    r <- tryCatch(httr::GET(paste0(url, path), query = query, httr::accept_json(), httr::timeout(60)),
                  error = function(e) NULL)
    if (!is.null(r) && httr::status_code(r) == 200) {
      return(jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8")))
    }
    if (!is.null(r) && httr::status_code(r) == 400) return(NULL)
    Sys.sleep(2 * i)
  }
  stop("Ensembl REST not reachable: ", path, call. = FALSE)
}

# LD proxies of each instrument, the instrument itself included
proxies <- map_dfr(unique(inst$rsid), function(id) {
  ld <- rest_get(sprintf("/ld/human/%s/%s", id, p$ld_population), list(r2 = p$r2, window_size = p$window_kb))
  Sys.sleep(0.1)
  other <- if (length(ld) && nrow(ld)) ld$variation2 else character()
  tibble(rsid = id, proxy = c(id, other), r2 = c(1, if (length(other)) as.numeric(ld$r2) else numeric()))
})

# VEP consequences of every proxy, 200 variants per request
vep_post <- function(ids) {
  for (i in 1:3) {
    r <- tryCatch(httr::POST(paste0(url, "/vep/human/id"), body = list(ids = ids), encode = "json",
                             httr::accept_json(), httr::timeout(300)), error = function(e) NULL)
    if (!is.null(r) && httr::status_code(r) == 200) {
      return(jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE))
    }
    Sys.sleep(5 * i)
  }
  stop("Ensembl VEP not reachable", call. = FALSE)
}
ids <- unique(proxies$proxy)
vep <- map_dfr(split(ids, ceiling(seq_along(ids) / 200)), function(chunk) {
  res <- vep_post(chunk)
  Sys.sleep(0.5)
  map_dfr(res, function(v) map_dfr(v$transcript_consequences %||% list(), function(tc) tibble(
    proxy = v$input %||% v$id,
    gene_symbol = tc$gene_symbol %||% NA_character_,
    consequence = paste(unlist(tc$consequence_terms), collapse = ","),
    amino_acids = tc$amino_acids %||% NA_character_,
    protein_start = as.integer(tc$protein_start %||% NA_integer_))))
})

# Protein-altering variants in the measured protein, the strongest-LD one per instrument
pav <- vep %>%
  filter(str_detect(consequence, paste(p$consequences, collapse = "|"))) %>%
  distinct(proxy, gene_symbol, consequence, amino_acids, protein_start)
flags <- inst %>%
  inner_join(proxies, by = "rsid", relationship = "many-to-many") %>%
  inner_join(pav, by = "proxy", relationship = "many-to-many") %>%
  filter(toupper(gene_symbol) == toupper(gene)) %>%
  arrange(desc(r2)) %>%
  slice_head(n = 1, by = c(gene, rsid)) %>%
  transmute(gene, rsid, pav_variant = proxy, pav_r2 = r2,
            pav_change = paste0(str_extract(consequence, paste(p$consequences, collapse = "|")),
                                if_else(is.na(amino_acids), "", paste0(" ", str_replace(amino_acids, "/", as.character(protein_start))))))

res <- inst %>%
  left_join(flags, by = c("gene", "rsid")) %>%
  mutate(protein_altering = !is.na(pav_variant)) %>%
  arrange(desc(protein_altering), gene)
write_tsv(res, file.path(out, "pqtl_protein_altering.tsv"))
print(filter(res, protein_altering), n = Inf)

write_counts("pqtl_protein_altering", c(
  proteins = n_distinct(res$gene),
  instruments = nrow(res),
  instruments_flagged = sum(res$protein_altering),
  proteins_flagged = n_distinct(res$gene[res$protein_altering])
))
