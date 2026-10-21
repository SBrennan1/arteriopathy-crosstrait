# Formats pairwise MTAG output for FUMA SNP2GENE (GRCh37 positions added).

source("R/setup.R")
out <- res_dir("mtag", "fuma_input")

ref <- read_univar("AA", cols = c("rsid", "chr", "pos", "effect_allele", "other_allele")) %>%
  distinct(rsid, .keep_all = TRUE)

for (t in names(cfg$mtag$outputs)) {
  m <- fread(file.path(cfg$mtag$dir, cfg$mtag$outputs[[t]]))
  m <- m %>%
    as_tibble() %>%
    transmute(rsid = SNP, effect_allele = toupper(A1), other_allele = toupper(A2),
              eaf = FRQ, beta = mtag_beta, se = mtag_se, pval = mtag_pval, N) %>%
    inner_join(ref %>% rename(ea_ref = effect_allele, oa_ref = other_allele), by = "rsid") %>%
    filter((effect_allele == ea_ref & other_allele == oa_ref) |
             (effect_allele == oa_ref & other_allele == ea_ref),
           grepl("^rs", rsid), !is.na(beta), !is.na(se), !is.na(pval)) %>%
    select(rsid, chr, pos, effect_allele, other_allele, eaf, beta, se, pval, N)

  message(t, ": ", nrow(m), " SNPs, ", sum(m$pval < 5e-8), " at P < 5e-8")
  fwrite(m, file.path(out, paste0(t, "_mtag.txt.gz")), sep = "\t")
}
