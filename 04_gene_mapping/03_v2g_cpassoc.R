# Variant-to-gene mapping of the CPASSOC lead SNPs. eQTL colocalization uses each trait
# with P < 0.05 at the lead (otherwise the most significant trait).

source("R/setup.R")
out <- res_dir("v2g")

loci <- read_tsv(file.path(cfg$results_dir, "cpassoc/cpassoc_loci_all.tsv"), show_col_types = FALSE)

p_lead <- as.matrix(loci[paste0("P_", traits)])
coloc_traits <- apply(p_lead, 1, function(p) {
  nominal <- which(p < 0.05)
  if (length(nominal) == 0) nominal <- which.min(p)
  paste(traits[nominal], collapse = ",")
})

leads <- loci %>% transmute(rsid, chr, pos, a1 = A1, a2 = A2, traits = coloc_traits)

res <- v2g_leads(leads, "CPASSOC")
v2g <- loci %>%
  select(rsid, P_SHet, cross_trait_tier) %>%
  inner_join(res$v2g, by = "rsid") %>%
  relocate(P_SHet, cross_trait_tier, .after = pos_hg38) %>%
  arrange(P_SHet)

write_tsv(v2g, file.path(out, "cpassoc_v2g.tsv"))
write_tsv(res$coloc, file.path(out, "cpassoc_gtex_coloc.tsv"))

cross <- filter(v2g, cross_trait_tier != "Single")
ev <- count(cross, evidence)
write_counts("v2g_cpassoc", c(
  lead_snps = nrow(v2g),
  cross_trait_loci = nrow(cross),
  setNames(ev$n, paste0("cross_trait_evidence_", str_replace_all(tolower(ev$evidence), "[^a-z0-9]+", "_"))),
  cross_trait_functional_support = sum(cross$evidence != "Nearest gene"),
  overridden = sum(v2g$overridden)
))
