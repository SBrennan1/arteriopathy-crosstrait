# Regional GTEx v10 eQTL subsets (gene body +/- 500 kb, GRCh38) from the all-pairs
# parquet files, for each tissue and gene with an FDR-significant eQTL MR result.

source("R/setup.R")
p <- cfg$qtl_coloc
out <- p$eqtl$dir
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- 500e3

existing <- list.files(out, pattern = "\\.tsv\\.gz$")
file_tissue <- function(t) p$eqtl$tissue_files[[t]] %||% t

hits <- read_tsv(file.path(cfg$results_dir, "mr/mr_eqtl.tsv"), show_col_types = FALSE) %>%
  filter(qval < p$fdr) %>%
  distinct(tissue, gene_id) %>%
  mutate(gtex = map_chr(tissue, file_tissue), gene = sub("\\..*", "", gene_id),
         have = map2_lgl(gtex, gene, ~ any(startsWith(existing, paste0(.x, ".", .y, ".")))))

# Gene bodies from the GENCODE annotation used by GTEx v10
genes <- fread(cmd = paste("gzip -dc", shQuote(cfg$gtex_v10$gencode), "| awk -F'\\t' '$3==\"gene\"'"),
               header = FALSE, sep = "\t", select = c(1, 4, 5, 9), col.names = c("chr", "start", "end", "attr"))
genes <- genes[!grepl("_PAR_Y", attr)][, .(gene = sub('.*gene_id "(ENSG\\d+).*', "\\1", attr),
                                           chr = sub("^chr", "", chr), start, end)]
win <- hits %>%
  left_join(as_tibble(genes), by = "gene") %>%
  mutate(from = pmax(1L, as.integer(start - w)), to = as.integer(end + w))

todo <- filter(win, !have)
message(nrow(hits), " tissue x gene pairs, ", nrow(todo), " without a subset")
if (any(is.na(todo$from))) {
  message("no GRCh38 gene body for ", sum(is.na(todo$from)), " genes; skipped")
  todo <- filter(todo, !is.na(from))
}

# One parquet read per tissue x chromosome
made <- todo %>%
  group_by(gtex, chr) %>%
  group_map(function(g, k) {
    f <- file.path(cfg$gtex_v10$dir, sprintf("%s.v10.allpairs.chr%s.parquet", k$gtex, k$chr))
    if (!file.exists(f)) {
      message("missing ", basename(f), " (", nrow(g), " genes)")
      return(tibble(gtex = k$gtex, gene = g$gene, n = NA_integer_))
    }
    d <- as.data.table(arrow::read_parquet(f))
    if ("phenotype_id" %in% names(d)) setnames(d, "phenotype_id", "gene_id")
    if ("start_distance" %in% names(d)) setnames(d, "start_distance", "tss_distance")
    d[, gene := sub("\\..*", "", gene_id)]
    d <- d[gene %chin% g$gene]
    d[, `:=`(chr = sub("_.*", "", variant_id), pos = as.integer(sub("^[^_]+_(\\d+)_.*", "\\1", variant_id)))]
    map_dfr(seq_len(nrow(g)), function(i) {
      r <- g[i, ]
      x <- d[gene == r$gene & pos >= r$from & pos <= r$to]
      keep <- intersect(c("chr", "pos", "gene_id", "variant_id", "tss_distance", "af", "ma_samples",
                          "ma_count", "pval_nominal", "slope", "slope_se"), names(x))
      if (nrow(x) > 0) {
        fwrite(x[, ..keep], file.path(out, sprintf("%s.%s.cis_chr%s_%d_%d.tsv.gz", k$gtex, r$gene, k$chr, r$from, r$to)),
               sep = "\t")
      }
      tibble(gtex = k$gtex, gene = r$gene, n = nrow(x))
    })
  }) %>%
  bind_rows()

print(made %>% summarise(genes = n(), written = sum(n > 0, na.rm = TRUE),
                         no_file = sum(is.na(n)), empty = sum(n == 0, na.rm = TRUE), .by = gtex))
write_counts("eqtl_subsets", c(pairs = nrow(hits), existing = sum(hits$have),
                               built = sum(made$n > 0, na.rm = TRUE), missing_parquet = sum(is.na(made$n))))
