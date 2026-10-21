# Genetic correlation between each arteriopathy and the matched negative controls. The
# null for a pair (A, B) is the rg of A with the controls matched to B, and of B with
# the controls matched to A. Run after 03_stream_ldsc.sh.

source("R/setup.R")
nc <- cfg$negative_controls
out <- res_dir("negative_controls")
mdir <- file.path(nc$dir, "munged")
dir.create(mdir, showWarnings = FALSE)

draws <- fread(nc$draws)
ld <- draws[ldsc_subset == TRUE]
ref <- if (file.exists(nc$reference_file)) fread(nc$reference_file) else data.table(label = character(), id = character(), h2_z_EUR = numeric())

ids <- unique(c(ld$id, ref$id))
have <- file.exists(file.path(nc$dir, "ldsc", paste0(ids, ".tsv.gz")))
if (any(!have)) warning(sum(!have), " LDSC controls not streamed yet, rerun 03_stream_ldsc.sh", call. = FALSE)
ids <- ids[have]
for (id in ids) munge(file.path(nc$dir, "ldsc", paste0(id, ".tsv.gz")), file.path(mdir, id))

# Every arteriopathy against every LDSC control
rg <- map_dfr(traits, function(t) {
  ldsc_rg(c(file.path(cfg$munged_dir, paste0(t, ".sumstats.gz")),
            file.path(mdir, paste0(ids, ".sumstats.gz"))),
          file.path(out, paste0("rg_", t, "_controls")))
}) %>%
  transmute(trait = p1, id = p2, rg, se, p, gcov_int, gcov_int_se)
rg_ref <- rg %>% inner_join(as_tibble(ref) %>% select(label, id), by = "id")
rg <- rg %>% filter(id %in% ld$id)
write_tsv(rg, file.path(out, "rg_controls.tsv"))

matched <- distinct(ld[id %in% ids, .(matched_to = trait, id)])
null_for <- function(a, b) {
  bind_rows(rg %>% filter(trait == a) %>% semi_join(filter(matched, matched_to == b), by = "id"),
            rg %>% filter(trait == b) %>% semi_join(filter(matched, matched_to == a), by = "id")) %>%
    filter(!is.na(rg))
}

# Pairs compared on |z| = |rg / se|; share of pairs with P < 0.05 compared with Fisher's
# exact test
obs <- read_tsv(file.path(cfg$results_dir, "ldsc/rg_pairwise.tsv"), show_col_types = FALSE)
emp_p <- function(null, x) (1 + sum(null >= x)) / (length(null) + 1)
pairs <- obs %>%
  transmute(trait1 = p1, trait2 = p2, rg, se, p, z = rg / se) %>%
  mutate(null = map2(trait1, trait2, null_for),
         n_null = map_int(null, nrow),
         null_median_abs_z = map_dbl(null, ~ median(abs(.x$rg / .x$se))),
         null_share_p05 = map_dbl(null, ~ mean(.x$p < 0.05)),
         empirical_p = map2_dbl(null, z, ~ emp_p(abs(.x$rg / .x$se), abs(.y))),
         null_median_abs_rg = map_dbl(null, ~ median(abs(.x$rg))),
         empirical_p_abs_rg = map2_dbl(null, rg, ~ emp_p(abs(.x$rg), abs(.y)))) %>%
  select(-null)
write_tsv(pairs, file.path(out, "rg_vs_controls.tsv"))
print(pairs %>% select(trait1, trait2, rg, z, n_null, null_share_p05, empirical_p))

ctrl_rg <- rg %>% filter(!is.na(p))
share <- fisher.test(matrix(c(sum(pairs$p < 0.05), sum(pairs$p >= 0.05),
                              sum(ctrl_rg$p < 0.05), sum(ctrl_rg$p >= 0.05)), 2))
rg_share <- tibble(arteriopathy_p05 = sum(pairs$p < 0.05), arteriopathy_pairs = nrow(pairs),
                   control_p05 = sum(ctrl_rg$p < 0.05), control_pairs = nrow(ctrl_rg),
                   odds_ratio = unname(share$estimate), fisher_p = share$p.value)
write_tsv(rg_share, file.path(out, "rg_share_significant.tsv"))
print(rg_share)

# Sample overlap check
message("Largest |gcov_int| between an arteriopathy and a control: ",
        round(max(abs(rg$gcov_int), na.rm = TRUE), 3))

rg_ref <- rg_ref %>%
  mutate(reference_h2_z = ref$h2_z_EUR[match(id, ref$id)]) %>%
  select(label, trait, id, rg, se, p, reference_h2_z)
write_tsv(rg_ref, file.path(out, "rg_reference.tsv"))
print(rg_ref)

write_counts("negative_controls_rg", c(
  ldsc_controls = length(setdiff(ids, ref$id)),
  arteriopathy_pairs_p05 = sum(pairs$p < 0.05),
  arteriopathy_pairs_beyond_null = sum(pairs$empirical_p < 0.05),
  control_pairs_p05 = sum(rg$p < 0.05, na.rm = TRUE),
  control_pairs = sum(!is.na(rg$p)),
  reference_pairs_p05 = sum(rg_ref$p < 0.05, na.rm = TRUE),
  rg_share_fisher_p = signif(rg_share$fisher_p, 3)
))
