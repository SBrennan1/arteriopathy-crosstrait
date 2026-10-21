# Summary counts for each stage, collected in results/counts.tsv
write_counts <- function(stage, x) {
  f <- file.path(cfg$results_dir, "counts.tsv")
  new <- tibble(stage = stage, item = names(x), n = unname(x))
  old <- if (file.exists(f)) read_tsv(f, show_col_types = FALSE) %>% filter(stage != !!stage) else NULL
  write_tsv(bind_rows(old, new), f)
  invisible(new)
}

# Harmonise trait labels to CeAD, FMD, IA, SCAD and AA
std_trait <- function(x) {
  key <- tolower(gsub("[^A-Za-z]", "", x))
  map <- c(cead = "CeAD", fmd = "FMD", fmdoutcome = "FMD", ia = "IA", scad = "SCAD",
           aa = "AA", ta = "AA", thoracicaneurysm = "AA", migraine = "Migraine")
  out <- unname(map[key])
  ifelse(is.na(out), x, out)
}
