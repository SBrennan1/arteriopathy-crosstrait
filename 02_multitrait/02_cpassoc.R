# CPASSOC (SHom and SHet; Zhu et al. 2015) on the same Z-score matrix as ASSET.
# Requires FunctionSet.R from the CPASSOC authors (config: cpassoc_functions).

source("R/setup.R")
source(cfg$cpassoc_functions)
select <- dplyr::select
out <- res_dir("cpassoc")

zm <- load_z_matrix(ref_trait = "CeAD")
corr <- read_gcov_int()

s_hom <- SHom(X = zm$z, SampleSize = zm$n, CorrMatrix = corr)
p_hom <- pchisq(s_hom, df = 1, lower.tail = FALSE)

# SHet null distribution estimated by simulation (seed set in setup.R)
para <- EstimateGamma(N = 1e4, SampleSize = zm$n, CorrMatrix = corr, correct = 1,
                      startCutoff = 0, endCutoff = 1, CutoffStep = 0.05,
                      isAllpossible = TRUE)
s_het <- SHet(X = zm$z, SampleSize = zm$n, CorrMatrix = corr)
p_het <- pgamma(s_het - para[3], shape = para[1], scale = para[2], lower.tail = FALSE)

cp <- tibble(rsid = rownames(zm$z), P_SHom = p_hom, P_SHet = p_het)
fwrite(cp, file.path(out, "cpassoc_all.tsv.gz"), sep = "\t")

leads <- clump_local(cp %>% select(rsid, pval = P_SHet))

loci <- cp %>%
  filter(rsid %in% leads$rsid) %>%
  left_join(zm$snps %>% transmute(rsid = SNP, A1, A2), by = "rsid") %>%
  left_join(snp_positions(leads$rsid), by = "rsid")

z_lead <- zm$z[loci$rsid, , drop = FALSE]
p_lead <- zm$p[loci$rsid, , drop = FALSE]

loci <- loci %>%
  bind_cols(nominal_support(z_lead, p_lead),
            as_tibble(z_lead) %>% rename_with(~ paste0("Z_", .x)),
            as_tibble(p_lead) %>% rename_with(~ paste0("P_", .x)),
            nearest_gene(loci$chr, loci$pos, load_genes()),
            min_trait_p = apply(p_lead, 1, min)) %>%
  mutate(
    # A: >= 2 traits nominal in the same direction
    # B: one trait genome-wide significant plus nominal support elsewhere (any direction)
    cross_trait_tier = case_when(
      n_pos >= 2 | n_neg >= 2 ~ "A",
      min_trait_p < 1e-8 & n_nominal >= 2 ~ "B",
      TRUE ~ "Single"
    )
  ) %>%
  add_hg38() %>%
  arrange(P_SHet)

write_tsv(loci, file.path(out, "cpassoc_loci_all.tsv"))
write_tsv(filter(loci, cross_trait_tier != "Single"), file.path(out, "cpassoc_loci.tsv"))

write_counts("cpassoc", c(
  snps_tested = nrow(cp),
  gws_shet = sum(cp$P_SHet < cfg$clump$p),
  independent_loci = nrow(loci),
  cross_trait_loci = sum(loci$cross_trait_tier != "Single")
))
