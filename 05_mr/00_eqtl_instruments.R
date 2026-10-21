# eQTL MR instruments: lead cis-eQTL per gene (GTEx v10, P < 5e-8), lifted to GRCh37.

source("R/setup.R")
p <- cfg$mr$eqtl
dir.create(p$instruments, recursive = TRUE, showWarnings = FALSE)
chain <- rtracklayer::import.chain(cfg$chain_hg38_to_hg19)

n <- imap_dfr(p$qcd, function(f, tissue) {
  d <- fread(f)
  setnames(d, "chr_col", "chr", skip_absent = TRUE)
  n_qcd <- nrow(d)
  d <- d[pval < p$p]
  d[, chr := sub("^chr", "", chr)]
  d[, pos := lift_eqtl_hg19(chr, pos, chain)]
  fwrite(d, file.path(p$instruments, paste0(tissue, ".tsv")), sep = "\t", na = "NA", quote = FALSE)
  tibble(tissue, qcd = n_qcd, instruments = nrow(d), genes = uniqueN(d$gene_id),
         no_pos37 = sum(is.na(d$pos)))
})
print(n)

write_counts("eqtl_instruments", set_names(n$instruments, n$tissue))
