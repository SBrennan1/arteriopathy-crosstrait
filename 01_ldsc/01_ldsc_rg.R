# Pairwise genetic correlation (LDSC) between the five arteriopathies, and the
# genetic covariance intercept matrix used by ASSET and CPASSOC.

source("R/setup.R")
out <- res_dir("ldsc")
dir.create(cfg$munged_dir, recursive = TRUE, showWarnings = FALSE)

for (t in traits) {
  munge(cfg$sumstats[[t]], file.path(cfg$munged_dir, t))
}
munged <- file.path(cfg$munged_dir, paste0(traits, ".sumstats.gz"))

rg <- map_dfr(seq_len(length(traits) - 1), function(i) {
  ldsc_rg(munged[i:length(traits)], file.path(out, paste0("rg_", traits[i])))
})
write_tsv(rg, file.path(out, "rg_pairwise.tsv"))

# Intercept matrix, checked against config/ldsc_gcov_int.tsv
gc <- diag(length(traits))
dimnames(gc) <- list(traits, traits)
for (k in seq_len(nrow(rg))) {
  gc[rg$p1[k], rg$p2[k]] <- gc[rg$p2[k], rg$p1[k]] <- rg$gcov_int[k]
}
write.table(round(gc, 4), file.path(out, "gcov_int.tsv"), sep = "\t", quote = FALSE,
            col.names = NA)

frozen <- read_gcov_int()
d <- max(abs(round(gc, 4) - frozen))
if (d > 1e-3) warning("gcov_int differs from config/ldsc_gcov_int.tsv (max diff ", d, ")")
