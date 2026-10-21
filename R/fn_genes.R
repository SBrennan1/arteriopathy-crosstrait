# Protein-coding genes from the Ensembl GRCh37 (release 87) GTF
load_genes <- function(path = cfg$gtf_grch37, protein_coding = TRUE) {
  g <- fread(cmd = paste("gzip -dc", shQuote(path), "| awk -F'\\t' '$3==\"gene\"'"),
             header = FALSE, sep = "\t",
             col.names = c("chr", "source", "feature", "start", "end", "score",
                           "strand", "frame", "attributes"))
  g %>%
    as_tibble() %>%
    mutate(
      chr = as.character(chr),
      gene_id = str_match(attributes, 'gene_id "([^"]+)"')[, 2],
      gene_name = str_match(attributes, 'gene_name "([^"]+)"')[, 2],
      biotype = str_match(attributes, 'gene_biotype "([^"]+)"')[, 2]
    ) %>%
    filter(!protein_coding | biotype == "protein_coding",
           !str_detect(gene_name, "^RP[0-9]|^AC[0-9]|^AL[0-9]|^LINC|^LOC[0-9]")) %>%
    select(chr, start, end, strand, gene_id, gene_name)
}

# Nearest gene body (distance 0 if the variant is inside the gene)
nearest_gene <- function(chr, pos, genes) {
  map2_dfr(as.character(chr), pos, function(ch, p) {
    g <- genes[genes$chr == ch, ]
    if (nrow(g) == 0) return(tibble(nearest_gene = NA_character_, nearest_distance = NA_integer_))
    d <- pmin(abs(p - g$start), abs(p - g$end))
    d[p >= g$start & p <= g$end] <- 0
    i <- which.min(d)
    tibble(nearest_gene = g$gene_name[i], nearest_distance = as.integer(d[i]))
  })
}
