# Analyses were run on GRCh37; output tables also carry GRCh38 coordinates.
add_hg38 <- function(d, chr_col = "chr", pos_col = "pos") {
  chain <- rtracklayer::import.chain(cfg$chain_hg19_to_hg38)
  key <- distinct(tibble(chr = as.character(d[[chr_col]]), pos = as.integer(d[[pos_col]])))
  gr <- GenomicRanges::GRanges(paste0("chr", key$chr), IRanges::IRanges(key$pos, width = 1))
  lo <- rtracklayer::liftOver(gr, chain)
  # keep only unique mappings
  ok <- lengths(lo) == 1
  key$chr_hg38 <- NA_character_
  key$pos_hg38 <- NA_integer_
  key$chr_hg38[ok] <- sub("^chr", "", as.character(GenomicRanges::seqnames(unlist(lo[ok]))))
  key$pos_hg38[ok] <- GenomicRanges::start(unlist(lo[ok]))

  d %>%
    mutate(.chr = as.character(.data[[chr_col]]), .pos = as.integer(.data[[pos_col]])) %>%
    left_join(key, by = c(.chr = "chr", .pos = "pos")) %>%
    select(-.chr, -.pos)
}
