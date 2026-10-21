# Formats the downloaded GWAS summary statistics (config: raw_sumstats) to a common
# layout (rsid chr pos effect_allele other_allele eaf beta se pval N, GRCh37).
# GRCh38 files are lifted to GRCh37; AAD odds ratios are converted to log(OR).

source("R/setup.R")
raw <- cfg$raw_sumstats
n_total <- function(t) cfg$traits[[t]]$ncase + cfg$traits[[t]]$ncontrol

# GRCh38 -> GRCh37. Variants without a unique mapping keep pos = NA and are matched
# by rsID downstream.
lift_to_hg19 <- function(d) {
  chain <- rtracklayer::import.chain(cfg$chain_hg38_to_hg19)
  lift_chr <- function(k) {
    lo <- rtracklayer::liftOver(
      GenomicRanges::GRanges(paste0("chr", k$chr), IRanges::IRanges(k$pos, width = 1)), chain)
    ok <- lengths(lo) == 1
    hit <- unlist(lo[ok])
    k$chr37 <- as.character(k$chr)
    k$pos37 <- NA_integer_
    k$chr37[ok] <- sub("^chr", "", as.character(GenomicRanges::seqnames(hit)))
    k$pos37[ok] <- GenomicRanges::start(hit)
    k
  }
  key <- d %>%
    distinct(chr, pos) %>%
    filter(!is.na(pos)) %>%
    group_split(chr) %>%
    map_dfr(lift_chr)
  d %>%
    left_join(key, by = c("chr", "pos")) %>%
    mutate(chr = coalesce(chr37, as.character(chr)), pos = pos37) %>%
    select(-chr37, -pos37)
}

# Column names as in each download; alleles read as text
alleles <- cols(effect_allele = "c", other_allele = "c")
read_raw <- list(
  CeAD = function(f) {
    vroom(f, show_col_types = FALSE, col_types = cols(Allele1 = "c", Allele2 = "c")) %>%
      rename(rsid = MarkerName, chr = CHR, pos = POS, effect_allele = Allele1,
             other_allele = Allele2, beta = Effect, se = StdErr, pval = pvalue) %>%
      mutate(N = n_total("CeAD"))
  },
  FMD = function(f) {
    vroom(f, show_col_types = FALSE, col_types = alleles) %>%
      rename(chr = chromosome, pos = base_pair_location, se = standard_error,
             eaf = effect_allele_frequency, pval = p_value) %>%
      mutate(N = n_cas + n_ctrls)
  },
  # N is the per-variant sample size given in the Bakker et al. file
  IA = function(f) {
    vroom(f, show_col_types = FALSE, col_types = cols(A_EFF = "c", A_NONEFF = "c")) %>%
      rename(chr = CHR, pos = BP, rsid = SNP, effect_allele = A_EFF, other_allele = A_NONEFF,
             eaf = Freq_EFF, beta = BETA, se = SE, pval = P)
  },
  SCAD = function(f) {
    vroom(f, show_col_types = FALSE, col_types = alleles) %>%
      rename(chr = chromosome, pos = base_pair_location, se = standard_error,
             eaf = effect_allele_frequency, pval = p_value, rsid = variant_id) %>%
      mutate(N = n_cases + n_ctrls)
  },
  # log(OR), with the SE recovered from the 95% CI
  AA = function(f) {
    vroom(f, show_col_types = FALSE, col_types = alleles) %>%
      rename(chr = chromosome, pos = base_pair_location, eaf = effect_allele_frequency,
             pval = p_value) %>%
      mutate(beta = log(odds_ratio),
             se = (log(ci_upper) - log(ci_lower)) / (2 * 1.96),
             se = if_else(is.finite(se), se, NA_real_),
             N = n_total("AA"))
  }
)

# Write ----
overview <- map_dfr(traits, function(t) {
  d <- read_raw[[t]](raw[[t]]$file)
  if (raw[[t]]$build == 38) d <- lift_to_hg19(d)
  # CeAD has no allele frequencies; palindromic SNPs are then dropped in MR harmonisation
  if (!"eaf" %in% names(d)) d$eaf <- NA_real_
  d <- d %>%
    mutate(effect_allele = toupper(effect_allele), other_allele = toupper(other_allele)) %>%
    select(rsid, chr, pos, effect_allele, other_allele, eaf, beta, se, pval, N)

  dir.create(dirname(cfg$sumstats[[t]]), recursive = TRUE, showWarnings = FALSE)
  fwrite(d, cfg$sumstats[[t]], sep = "\t", na = "NA", quote = FALSE)

  tibble(trait = t, variants = nrow(d), no_position = sum(is.na(d$pos)),
         N_min = min(d$N, na.rm = TRUE), N_max = max(d$N, na.rm = TRUE))
})
write_tsv(overview, file.path(res_dir("sumstats"), "sumstats_summary.tsv"))

write_counts("sumstats", c(
  set_names(overview$variants, paste0(overview$trait, "_variants")),
  set_names(overview$N_max, paste0(overview$trait, "_N_max")),
  AA_no_position = overview$no_position[overview$trait == "AA"]
))
