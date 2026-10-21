# Genome-wide munged summary statistics for ASSET, CPASSOC and MTAG, restricted to
# SNPs in the 1000 Genomes EUR reference. Allele frequency and position are added for
# MTAG (reference frequencies for CeAD).

source("R/setup.R")
dir.create(cfg$munged_full_dir, recursive = TRUE, showWarnings = FALSE)

# Value of `col` in d for each row of m, matching rsID and alleles in either orientation
lookup <- function(m, d, col, flip = identity) {
  d <- d[rsid %chin% m$SNP]
  key <- unique(rbind(
    d[, .(SNP = rsid, A1 = effect_allele, A2 = other_allele, v = get(col))],
    d[, .(SNP = rsid, A1 = other_allele, A2 = effect_allele, v = flip(get(col)))]
  ), by = c("SNP", "A1", "A2"))
  key[m, on = c("SNP", "A1", "A2"), v]
}

ref_frq <- function() {
  if (!file.exists(paste0(cfg$ld_ref, ".bed"))) stop("1000G EUR panel not found at ", cfg$ld_ref, ".bed")
  f <- file.path(tempdir(), "ref_frq")
  system2(plink_bin(), c("--bfile", cfg$ld_ref, "--freq", "--out", f), stdout = FALSE)
  fread(paste0(f, ".frq"))[, .(rsid = SNP, effect_allele = A1, other_allele = A2, eaf = MAF)]
}

panel <- fread(paste0(cfg$ld_ref, ".bim"), select = 2)$V2

n <- map_dfr(traits, function(t) {
  out_file <- file.path(cfg$munged_full_dir, paste0(t, ".sumstats.gz"))
  if (!file.exists(out_file)) {
    tmp <- file.path(cfg$munged_full_dir, paste0(t, "_munge"))
    munge(cfg$sumstats[[t]], tmp, hm3 = FALSE)
    m <- fread(paste0(tmp, ".sumstats.gz"))[!is.na(Z) & SNP %chin% panel]

    u <- as.data.table(read_univar(t, cols = c("rsid", "pos", "effect_allele", "other_allele", "eaf")))
    frq <- if (all(is.na(u$eaf))) ref_frq() else u
    m[, FRQ := lookup(m, frq, "eaf", flip = function(x) 1 - x)]
    m[, BP := lookup(m, u, "pos")]
    # MTAG's munge step requires a P column
    m[, P := 2 * pnorm(-abs(Z))]

    fwrite(m, out_file, sep = "\t", na = "NA")
    file.remove(paste0(tmp, ".sumstats.gz"))
  }
  m <- fread(out_file, select = c("FRQ", "BP"))
  tibble(trait = t, snps = nrow(m), no_frq = sum(is.na(m$FRQ)), no_bp = sum(is.na(m$BP)))
})
print(n)

write_counts("munge_full", c(set_names(n$snps, paste0(n$trait, "_snps")),
                             set_names(n$no_frq, paste0(n$trait, "_no_frq"))))
