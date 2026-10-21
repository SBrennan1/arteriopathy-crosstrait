# Cross-trait colocalization at the ASSET loci: coloc.abf (primary) with coloc.susie
# where both traits have credible sets.

source("R/setup.R")
out <- res_dir("coloc")

pp <- c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
abf <- read_tsv(file.path(out, "asset_coloc_abf.tsv"), show_col_types = FALSE) %>%
  select(locus_rsid, trait1, trait2, n_snps, all_of(pp))
susie <- read_tsv(file.path(out, "asset_susie_coloc.tsv"), show_col_types = FALSE) %>%
  transmute(locus_rsid, trait1, trait2, susie_PP.H3 = PP.H3.abf, susie_PP.H4 = PP.H4.abf, hit1, hit2)

coloc <- abf %>%
  left_join(susie, by = c("locus_rsid", "trait1", "trait2")) %>%
  left_join(read_tsv(file.path(cfg$results_dir, "asset/asset_loci.tsv"), show_col_types = FALSE) %>%
              select(locus_rsid = rsid, chr, pos, chr_hg38, pos_hg38, nearest_gene),
            by = "locus_rsid") %>%
  mutate(method = "coloc.abf",
         class = case_when(
           is.na(PP.H4.abf) ~ "Not tested",
           PP.H4.abf >= 0.8 ~ "Strong",
           PP.H4.abf >= 0.5 ~ "Moderate",
           PP.H3.abf > 0.5 ~ "Distinct",
           TRUE ~ "None"
         ),
         susie_class = case_when(
           is.na(susie_PP.H4) ~ NA_character_,
           susie_PP.H4 >= 0.8 ~ "Strong",
           susie_PP.H4 >= 0.5 ~ "Moderate",
           susie_PP.H3 > 0.5 ~ "Distinct",
           TRUE ~ "None"
         )) %>%
  arrange(chr, pos, trait1, trait2)

extra <- anti_join(susie, abf, by = c("locus_rsid", "trait1", "trait2"))
if (nrow(extra)) message(nrow(extra), " SuSiE pairs without a coloc.abf result")

write_tsv(coloc, file.path(out, "cross_trait_coloc.tsv"))

strong <- filter(coloc, class == "Strong")
write_counts("cross_trait_coloc", c(
  pairs_tested = sum(!is.na(coloc$PP.H4.abf)),
  pairs_h4_0.8 = nrow(strong),
  pairs_h4_0.5_0.8 = sum(coloc$class == "Moderate"),
  pairs_distinct = sum(coloc$class == "Distinct"),
  loci_h4_0.8 = n_distinct(strong$locus_rsid),
  pairs_susie = sum(!is.na(coloc$susie_PP.H4)),
  pairs_susie_h4_0.8 = sum(coloc$susie_class == "Strong", na.rm = TRUE)
))
print(count(coloc, class))
print(coloc %>% filter(!is.na(susie_PP.H4)) %>%
        transmute(nearest_gene, trait1, trait2, abf_H4 = round(PP.H4.abf, 2), susie_H4 = round(susie_PP.H4, 2),
                  susie_H3 = round(susie_PP.H3, 2)), n = Inf)
