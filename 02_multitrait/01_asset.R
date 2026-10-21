# ASSET subset-based meta-analysis across the five arteriopathies, with sample overlap
# handled through the LDSC genetic covariance intercepts.

source("R/setup.R")
library(ASSET)
out <- res_dir("asset")

zm <- load_z_matrix(ref_trait = "CeAD")
message(nrow(zm$z), " SNPs common to all traits")

# ASSET on Z-scores (beta = Z, se = 1), restricted to SNPs with P < 0.05 in any trait
keep <- which(rowSums(zm$p < 0.05) > 0)

fit <- h.traits(
  snp.vars = rownames(zm$z)[keep],
  traits.lab = traits,
  beta.hat = zm$z[keep, ],
  sigma.hat = matrix(1, length(keep), length(traits), dimnames = list(NULL, traits)),
  cor = read_gcov_int(),
  ncase = map_dbl(cfg$traits, "ncase"),
  ncntl = map_dbl(cfg$traits, "ncontrol"),
  cor.numr = FALSE,
  side = 2,
  search = 2
)

s2 <- fit$Subset.2sided
subset_label <- function(m) apply(m, 1, function(x) paste(colnames(m)[x], collapse = ","))

asset <- tibble(
  rsid = fit$snp.vars,
  asset_p = s2$pval,
  asset_beta_pos = s2$beta.1,
  asset_beta_neg = s2$beta.2,
  asset_subset_pos = subset_label(s2$pheno.1),
  asset_subset_neg = subset_label(s2$pheno.2)
)
fwrite(asset, file.path(out, "asset_all.tsv.gz"), sep = "\t")

# Independent genome-wide significant loci ----
leads <- clump_local(asset %>% select(rsid, pval = asset_p))

genes <- load_genes()

loci <- asset %>%
  filter(rsid %in% leads$rsid) %>%
  left_join(zm$snps %>% transmute(rsid = SNP, A1, A2), by = "rsid") %>%
  left_join(snp_positions(leads$rsid), by = "rsid")

z_lead <- zm$z[loci$rsid, , drop = FALSE]
p_lead <- zm$p[loci$rsid, , drop = FALSE]

loci <- loci %>%
  bind_cols(nominal_support(z_lead, p_lead),
            as_tibble(z_lead) %>% rename_with(~ paste0("Z_", .x)),
            as_tibble(p_lead) %>% rename_with(~ paste0("P_", .x)),
            nearest_gene(loci$chr, loci$pos, genes)) %>%
  add_hg38() %>%
  arrange(asset_p)

# Pleiotropic loci: nominal support (P < 0.05) in at least two traits
pleio <- loci %>% filter(n_nominal >= 2)

write_tsv(loci, file.path(out, "asset_loci_all.tsv"))
write_tsv(pleio, file.path(out, "asset_loci.tsv"))

write_counts("asset", c(
  common_snps = nrow(zm$z),
  snps_tested = length(keep),
  gws_snps = sum(asset$asset_p < cfg$clump$p),
  independent_loci = nrow(loci),
  pleiotropic_loci = nrow(pleio),
  bidirectional = sum(pleio$effect_pattern == "Bidirectional")
))
