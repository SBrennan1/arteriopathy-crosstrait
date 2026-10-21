# Overlap between ASSET, CPASSOC and pairwise MTAG loci. Leads from different methods
# are merged into one locus if within 500 kb or in LD (r2 > 0.1, 1000G EUR).

source("R/setup.R")
out <- res_dir("overlap")

leads <- bind_rows(
  read_tsv(file.path(cfg$results_dir, "asset/asset_loci.tsv"), show_col_types = FALSE) %>%
    transmute(method = "ASSET", rsid, chr, pos, p = asset_p),
  read_tsv(file.path(cfg$results_dir, "cpassoc/cpassoc_loci.tsv"), show_col_types = FALSE) %>%
    transmute(method = "CPASSOC", rsid, chr, pos, p = P_SHet),
  read_tsv(file.path(cfg$results_dir, "mtag/mtag_leads.tsv"), show_col_types = FALSE) %>%
    transmute(method = "MTAG", rsid, chr, pos, p = p_mtag)
) %>%
  group_by(method, rsid) %>% slice_min(p, n = 1, with_ties = FALSE) %>% ungroup()

# Single-linkage grouping on distance or LD
group_chr <- function(d) {
  link <- abs(outer(d$pos, d$pos, "-")) <= 500e3
  r2 <- ld_r2(d$rsid)
  if (!is.null(r2)) {
    i <- match(d$rsid, rownames(r2))
    have <- which(!is.na(i))
    link[have, have] <- link[have, have] | r2[i[have], i[have]] > 0.1
  }
  diag(link) <- FALSE
  g <- igraph::graph_from_adjacency_matrix(link * 1, mode = "undirected")
  igraph::components(g)$membership
}

leads <- leads %>%
  arrange(chr, pos) %>%
  group_by(chr) %>%
  mutate(locus = paste0(chr, "_", group_chr(pick(rsid, pos)))) %>%
  ungroup()

loci <- leads %>%
  group_by(locus) %>%
  summarise(chr = first(chr), start = min(pos), end = max(pos),
            lead_rsid = rsid[which.min(p)],
            ASSET = any(method == "ASSET"), CPASSOC = any(method == "CPASSOC"),
            MTAG = any(method == "MTAG"),
            rsids = paste(unique(rsid), collapse = ","), .groups = "drop") %>%
  mutate(methods = pmap_chr(list(ASSET, CPASSOC, MTAG), function(a, c, m)
    paste(c("ASSET", "CPASSOC", "MTAG")[c(a, c, m)], collapse = "+")))

# Univariate P at each locus lead
univ <- map_dfr(traits, function(t) {
  read_univar(t, cols = c("rsid", "pval")) %>% filter(rsid %in% loci$lead_rsid) %>%
    mutate(trait = t)
}) %>%
  pivot_wider(names_from = trait, values_from = pval, names_prefix = "P_")

loci <- loci %>% left_join(univ, by = c("lead_rsid" = "rsid"))

summary_tab <- loci %>% count(methods, name = "n_loci") %>% arrange(desc(n_loci))

pairwise <- combn(c("ASSET", "CPASSOC", "MTAG"), 2, simplify = FALSE) %>%
  map_dfr(function(m) {
    a <- loci[[m[1]]]; b <- loci[[m[2]]]
    tibble(method1 = m[1], method2 = m[2], n1 = sum(a), n2 = sum(b), shared = sum(a & b),
           pct_of_1 = round(100 * shared / n1, 1), pct_of_2 = round(100 * shared / n2, 1),
           jaccard = round(shared / sum(a | b), 3))
  })

write_tsv(loci, file.path(out, "locus_by_method.tsv"))
write_tsv(summary_tab, file.path(out, "method_combinations.tsv"))
write_tsv(pairwise, file.path(out, "method_pairwise_overlap.tsv"))

write_counts("overlap", c(
  total_loci = nrow(loci),
  all_three = sum(loci$ASSET & loci$CPASSOC & loci$MTAG),
  setNames(pairwise$shared, paste0(pairwise$method1, "_", pairwise$method2))
))
