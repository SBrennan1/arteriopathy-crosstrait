# Cytogenetic band of each ASSET lead variant (Ensembl GRCh37 REST).

source("R/setup.R")

loci <- read_stage("asset", "asset_loci.tsv") %>%
  select(rsid, chr, pos, chr_hg38, pos_hg38, nearest_gene)

loci %>%
  mutate(cytoband = map2_chr(chr, pos, ensembl_cytoband)) %>%
  write_tsv(file.path(res_dir("asset"), "asset_cytobands.tsv"))
