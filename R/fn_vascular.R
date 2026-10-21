# Helpers for 11_vascular_expression/

read_expr <- function(dataset, min_cells = cfg$vascular$min_cells) {
  read_stage("vascular", paste0("expr_", dataset, ".tsv.gz")) %>%
    mutate(across(c(donor, region, cell_type), as.character)) %>%
    filter(n_cells >= min_cells, !cell_type %in% c("unknown", "Unknown", "NA"))
}

# Donor-level means, then the average over donors. by: the grouping kept besides gene
donor_average <- function(e, by) {
  e %>%
    summarise(mean_cp10k = weighted.mean(mean_cp10k, n_cells), mean_log = weighted.mean(mean_log, n_cells),
              pct = weighted.mean(pct, n_cells), n_cells = sum(n_cells),
              .by = all_of(c("donor", by, "gene_id", "gene"))) %>%
    summarise(mean_cp10k = mean(mean_cp10k), mean_log = mean(mean_log), pct = mean(pct),
              n_donors = n(), n_cells = sum(n_cells), .by = all_of(c(by, "gene_id", "gene")))
}

# Specificity of each gene for each cell type: its mean expression in that cell type
# over the sum across cell types; tau summarises how restricted the gene is overall
celltype_profile <- function(e) {
  donor_average(e, "cell_type") %>%
    mutate(specificity = coalesce(mean_cp10k / sum(mean_cp10k), 0),
           tau = tau_index(mean_log),
           .by = gene_id)
}

tau_index <- function(x) {
  if (length(x) < 2 || max(x) <= 0) return(rep(NA_real_, length(x)))
  rep(sum(1 - x / max(x)) / (length(x) - 1), length(x))
}

# One unit per target gene: the target and the other genes of every locus it was
# chosen at. Genes that are a target anywhere in the set are left out of the pools.
make_units <- function(lg) {
  targets <- unique(lg$gene[lg$is_target])
  tl <- lg %>% filter(is_target) %>% distinct(locus, source, unit = gene)
  lg %>%
    inner_join(tl, by = c("locus", "source"), relationship = "many-to-many") %>%
    filter(gene == unit | !gene %in% targets) %>%
    transmute(unit, gene, is_target = gene == unit) %>%
    distinct()
}

# Locus-matched permutation: is the mean score of the target genes higher than that of
# one gene drawn at random from each of the same loci? d has unit, is_target, score.
locus_perm <- function(d, n_perm = cfg$vascular$n_perm) {
  d <- d %>%
    filter(!is.na(score)) %>%
    filter(any(is_target), any(!is_target), .by = unit)
  if (nrow(d) == 0) return(tibble(n_units = 0L))
  tgt <- d %>% filter(is_target) %>% distinct(unit, .keep_all = TRUE)
  other <- d %>% filter(!is_target)
  pools <- split(other$score, other$unit)[tgt$unit]
  
  null <- vapply(pools, function(x) x[sample.int(length(x), n_perm, replace = TRUE)], numeric(n_perm))
  nm <- rowMeans(matrix(null, nrow = n_perm))
  obs <- mean(tgt$score)
  pctl <- map2_dbl(tgt$score, pools, ~ (sum(.y < .x) + 0.5 * sum(.y == .x)) / length(.y))
  tibble(n_units = nrow(tgt), n_other = nrow(other), observed = obs,
         null_mean = mean(nm), null_lo = unname(quantile(nm, 0.025)),
         null_hi = unname(quantile(nm, 0.975)),
         p_greater = (1 + sum(nm >= obs)) / (1 + n_perm),
         p_less = (1 + sum(nm <= obs)) / (1 + n_perm),
         mean_percentile = mean(pctl))
}

# Expression-matched permutation (sensitivity to locus_perm): each target gene is
# compared with non-target locus genes, from any locus, in the same decile of overall
# expression, so that a target set of highly expressed genes is not credited with
# specificity that comes from expression level alone. d has gene, is_target, score,
# level (overall expression).
expr_matched_perm <- function(d, n_perm = cfg$vascular$n_perm, n_bins = 10) {
  d <- d %>% filter(!is.na(score), !is.na(level)) %>% distinct(gene, .keep_all = TRUE)
  d$bin <- cut(rank(d$level, ties.method = "first"), n_bins, labels = FALSE)
  tgt <- d %>% filter(is_target)
  other <- d %>% filter(!is_target)
  pools <- split(other$score, other$bin)[as.character(tgt$bin)]
  ok <- lengths(pools) > 0
  tgt <- tgt[ok, ]
  pools <- pools[ok]
  if (nrow(tgt) == 0) return(tibble(n_units = 0L))
  null <- vapply(pools, function(x) x[sample.int(length(x), n_perm, replace = TRUE)], numeric(n_perm))
  nm <- rowMeans(matrix(null, nrow = n_perm))
  obs <- mean(tgt$score)
  tibble(n_units = nrow(tgt), n_other = nrow(other), observed = obs,
         null_mean = mean(nm), null_lo = unname(quantile(nm, 0.025)),
         null_hi = unname(quantile(nm, 0.975)),
         p_greater = (1 + sum(nm >= obs)) / (1 + n_perm),
         p_less = (1 + sum(nm <= obs)) / (1 + n_perm))
}

# Benjamini-Hochberg q-values for both tails within each test family
add_q <- function(d, ...) {
  d %>% mutate(q_greater = p.adjust(p_greater, "BH"), q_less = p.adjust(p_less, "BH"), .by = c(...))
}

# Locus genes with every locus re-assigned to its nearest gene (an assignment that uses
# no expression data), for the sensitivity analyses
nearest_gene_targets <- function(lg, loci) {
  lg %>%
    left_join(loci %>% distinct(locus, source, nearest_gene), by = c("locus", "source")) %>%
    mutate(is_target = coalesce(gene == nearest_gene, FALSE)) %>%
    select(-nearest_gene)
}

# Units for loci from the given sources, optionally only loci involving a trait
units_for <- function(lg, loci, sources, trait = NULL) {
  keep <- loci %>% filter(source %in% sources)
  if (!is.null(trait)) keep <- keep %>% filter(map_lgl(str_split(traits, ","), ~ trait %in% .x))
  lg %>% semi_join(keep, by = c("locus", "source")) %>% make_units()
}
