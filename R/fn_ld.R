plink_bin <- function() genetics.binaRies::get_plink_binary()

# PLINK clumping against the local 1000G EUR panel. d needs rsid and pval.
clump_local <- function(d, p = cfg$clump$p, r2 = cfg$clump$r2, kb = cfg$clump$kb) {
  d <- d %>% filter(pval < p) %>% select(rsid, pval)
  if (nrow(d) == 0) return(d)
  ieugwasr::ld_clump(d, clump_kb = kb, clump_r2 = r2, clump_p = p,
                     bfile = cfg$ld_ref, plink_bin = plink_bin())
}

# Pairwise r2 among a set of rsIDs (SNPs absent from the panel are dropped)
ld_r2 <- function(rsids, bfile = cfg$ld_ref) {
  rsids <- unique(rsids)
  if (length(rsids) < 2) return(NULL)
  r <- tryCatch(ieugwasr::ld_matrix(rsids, with_alleles = FALSE, bfile = bfile,
                                    plink_bin = plink_bin()),
                error = function(e) NULL)
  if (is.null(r)) NULL else r^2
}
