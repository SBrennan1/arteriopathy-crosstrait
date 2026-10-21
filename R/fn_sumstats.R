read_univar <- function(trait, cols = NULL) {
  d <- vroom(cfg$sumstats[[trait]], show_col_types = FALSE, col_select = all_of(cols %||% everything()))
  d %>%
    mutate(across(any_of(c("effect_allele", "other_allele")), toupper),
           across(any_of(c("chr", "pos")), as.integer))
}

read_munged <- function(trait) {
  d <- fread(file.path(cfg$munged_full_dir, paste0(trait, ".sumstats.gz")),
             select = c("SNP", "A1", "A2", "Z", "N"))
  d <- d[!is.na(Z) & !duplicated(SNP)]
  d[, P := 2 * pnorm(-abs(Z))]
  d
}

# Put Z on the reference allele. Drops SNPs whose alleles don't match the
# reference in either orientation rather than silently flipping them.
align_z <- function(d, ref) {
  d <- merge(d, ref, by = "SNP")
  same <- d$A1 == d$A1_ref & d$A2 == d$A2_ref
  swap <- d$A1 == d$A2_ref & d$A2 == d$A1_ref
  if (any(!(same | swap))) message("  dropping ", sum(!(same | swap)), " SNPs with mismatched alleles")
  d <- d[same | swap]
  d[!(A1 == A1_ref), Z := -Z]
  d[, `:=`(A1 = A1_ref, A2 = A2_ref, A1_ref = NULL, A2_ref = NULL)]
  d[, P := 2 * pnorm(-abs(Z))]
  d
}

# Load all traits from munged files, align to one reference, keep SNPs present
# in every trait. Returns list(snps, z, p) with z/p as SNP x trait matrices.
load_z_matrix <- function(ref_trait = "CeAD") {
  ss <- setNames(lapply(traits, read_munged), traits)
  ref <- ss[[ref_trait]][, .(SNP, A1_ref = A1, A2_ref = A2)]
  ss <- lapply(ss, align_z, ref = ref)
  common <- sort(Reduce(intersect, lapply(ss, `[[`, "SNP")))
  ss <- lapply(ss, function(d) d[match(common, SNP)])
  stopifnot(all(vapply(ss, function(d) identical(d$SNP, common), logical(1))))

  z <- sapply(ss, `[[`, "Z")
  rownames(z) <- common
  list(snps = ss[[ref_trait]][, .(SNP, A1, A2)], z = z, p = 2 * pnorm(-abs(z)),
       n = sapply(ss, function(d) max(d$N)))
}

read_gcov_int <- function() {
  m <- read.delim(cfg$gcov_int, row.names = 1, check.names = FALSE)
  m <- as.matrix(m)[traits, traits]
  stopifnot(isSymmetric(m))
  m
}

# GRCh37 positions for a set of rsIDs. AA first (best coverage), then the others
snp_positions <- function(rsids) {
  found <- tibble(rsid = character(), chr = integer(), pos = integer())
  for (t in c("AA", setdiff(traits, "AA"))) {
    todo <- setdiff(rsids, found$rsid)
    if (length(todo) == 0) break
    found <- bind_rows(found, read_univar(t, cols = c("rsid", "chr", "pos")) %>%
                         filter(rsid %in% todo) %>% distinct(rsid, .keep_all = TRUE))
  }
  found
}

# Which traits support a variant at P < alpha, split by direction of Z.
# z and p are data frames / matrices with one column per trait.
nominal_support <- function(z, p, alpha = 0.05) {
  z <- as.matrix(z); p <- as.matrix(p)
  sig <- p < alpha & !is.na(p)
  pos <- sig & z > 0
  neg <- sig & z < 0
  collapse <- function(m) apply(m, 1, function(x) {
    s <- paste(colnames(m)[x], collapse = ",")
    if (s == "") NA_character_ else s
  })
  tibble(
    n_nominal = rowSums(sig),
    n_pos = rowSums(pos),
    n_neg = rowSums(neg),
    traits_pos = collapse(pos),
    traits_neg = collapse(neg),
    effect_pattern = case_when(
      n_pos >= 1 & n_neg >= 1 ~ "Bidirectional",
      n_pos >= 2 ~ "Concordant positive",
      n_neg >= 2 ~ "Concordant negative",
      TRUE ~ "Single trait"
    )
  )
}
