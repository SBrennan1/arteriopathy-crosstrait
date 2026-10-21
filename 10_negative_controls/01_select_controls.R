# Negative controls: Pan-UKB phecodes matched to each arteriopathy on EUR case count,
# drawn with a fixed seed (saved to config/ once drawn). Also writes the ASSET regions
# and HapMap3 position map used by 02_extract_regions.sh and 03_stream_ldsc.sh.

source("R/setup.R")
nc <- cfg$negative_controls
dir.create(nc$dir, recursive = TRUE, showWarnings = FALSE)
out <- res_dir("negative_controls")

# Pan-UKB phenotype manifest
read_manifest <- function() {
  f <- file.path(nc$dir, "phenotype_manifest.tsv.bgz")
  if (!file.exists(f)) download.file(nc$manifest, f, mode = "wb")
  man <- fread(cmd = paste("gzip -dc", shQuote(f)), na.strings = c("", "NA"),
               colClasses = list(character = "phenocode"))
  if (!all(c("aws_link", "aws_link_tabix") %in% names(man))) {
    s3_col <- function(pattern) {
      hit <- names(man)[map_lgl(man, ~ is.character(.x) && any(grepl(pattern, head(.x, 100))))]
      if (length(hit) == 0) stop("no manifest column matching ", pattern)
      hit[1]
    }
    to_https <- function(x) sub("^s3://([^/]+)/", "https://\\1.s3.amazonaws.com/", x)
    man[, aws_link := to_https(get(s3_col("^s3://.*\\.tsv\\.bgz$")))]
    man[, aws_link_tabix := to_https(get(s3_col("^s3://.*\\.tbi$")))]
  }
  man[, id := sub("\\.tsv\\.bgz$", "", basename(aws_link))]
  man
}

# Draws ----
if (!file.exists(nc$draws)) {
  man <- read_manifest()
  need <- c("trait_type", "phenocode", "pheno_sex", "description", "category", "phenotype_qc_EUR",
            "n_cases_EUR", "n_controls_EUR", "sldsc_25bin_h2_z_EUR", "aws_link", "aws_link_tabix")
  miss <- setdiff(need, names(man))
  if (length(miss)) stop("manifest is missing columns: ", paste(miss, collapse = ", "))

  pool <- man[trait_type == "phecode" & pheno_sex == "both_sexes" & phenotype_qc_EUR %in% nc$qc_keep]
  drop <- grepl(paste(nc$exclude, collapse = "|"), paste(pool$category, pool$description),
                ignore.case = TRUE)
  message(nrow(pool), " phecodes pass EUR QC, ", sum(!drop), " after exclusions")
  fwrite(pool[drop, .(phenocode, description, category, n_cases_EUR)][order(category, description)],
         file.path(out, "excluded_phecodes.tsv"), sep = "\t")
  pool <- pool[!drop]
  tol <- function(t) if (is.list(nc$case_tol)) nc$case_tol[[t]] %||% nc$case_tol$default else nc$case_tol
  eligible <- rbindlist(lapply(traits, function(t) {
    n_case <- cfg$traits[[t]]$ncase
    pool[abs(n_cases_EUR - n_case) <= tol(t) * n_case,
         .(trait = t, phenocode, description, category, n_cases_EUR)]
  }))
  fwrite(eligible[order(trait, category, description)], file.path(out, "eligible_controls.tsv"), sep = "\t")
  stopifnot(!anyDuplicated(pool$id))

  if (nc$match_h2_z) {
    h2 <- read_tsv(file.path(cfg$results_dir, "ldsc/h2_liability.tsv"), show_col_types = FALSE)
    h2_z <- setNames(h2$h2_liab / h2$h2_liab_se, h2$trait)
  }

  draws <- rbindlist(lapply(seq_along(traits), function(i) {
    t <- traits[i]
    n_case <- cfg$traits[[t]]$ncase
    cand <- pool[abs(n_cases_EUR - n_case) <= tol(t) * n_case]
    if (nc$match_h2_z) {
      cand <- cand[abs(sldsc_25bin_h2_z_EUR - h2_z[[t]]) <= nc$h2z_tol * abs(h2_z[[t]])]
    }
    n <- min(nc$n_draws, nrow(cand))
    if (n < nc$n_draws) warning(t, ": only ", n, " matched controls", call. = FALSE)
    if (n == 0) return(NULL)
    set.seed(nc$seed + i)
    cand[sample(.N, n)][, `:=`(trait = t, draw = seq_len(.N), ldsc_subset = seq_len(.N) <= nc$n_ldsc)]
  }))

  draws <- draws[, .(trait, draw, ldsc_subset, id, phenocode, description, category,
                     n_cases_EUR, n_controls_EUR, h2_z_EUR = sldsc_25bin_h2_z_EUR,
                     aws_link, aws_link_tabix)]
  fwrite(draws, nc$draws, sep = "\t")
}
draws <- fread(nc$draws)

# Reference trait (hair color) ----
if (!file.exists(nc$reference_file)) {
  man <- read_manifest()
  ref <- imap_dfr(nc$reference, function(r, label) {
    hit <- man[trait_type == r$trait_type & phenocode == as.character(r$phenocode) & pheno_sex == "both_sexes"]
    if (!is.null(r$coding)) hit <- hit[as.character(coding) == as.character(r$coding)]
    if (nrow(hit) != 1) {
      print(man[phenocode == as.character(r$phenocode),
                .SD, .SDcols = intersect(c("trait_type", "pheno_sex", "coding", "coding_description",
                                           "description", "n_cases_EUR"), names(man))])
      stop(label, ": ", nrow(hit), " manifest rows match; set trait_type and coding in config.yml")
    }
    what <- if ("coding_description" %in% names(hit)) paste(hit$description, hit$coding_description, sep = ": ") else hit$description
    hit[, .(label, trait = "reference", draw = NA_integer_, ldsc_subset = TRUE, id, phenocode,
            description = what, category, n_cases_EUR, n_controls_EUR,
            h2_z_EUR = sldsc_25bin_h2_z_EUR, aws_link, aws_link_tabix)]
  })
  fwrite(ref, nc$reference_file, sep = "\t")
}
ref <- fread(nc$reference_file)
print(ref[, .(label, id, description, n_cases_EUR, n_controls_EUR, h2_z_EUR)])

smry <- draws[, .(n_drawn = .N, n_ldsc = sum(ldsc_subset), median_cases = median(n_cases_EUR),
                 min_cases = min(n_cases_EUR), max_cases = max(n_cases_EUR),
                 median_h2_z = median(h2_z_EUR, na.rm = TRUE)), by = trait]
smry[, arteriopathy_cases := map_dbl(trait, ~ cfg$traits[[.x]]$ncase)]
print(smry)
message(uniqueN(draws$id), " distinct controls (", uniqueN(draws[ldsc_subset == TRUE]$id),
        " in the LDSC subset)")
fwrite(smry, file.path(out, "draws_summary.tsv"), sep = "\t")

# Regions for 02_extract_regions.sh ----
w <- cfg$coloc$window_kb * 1e3
loci <- read_tsv(file.path(cfg$results_dir, "asset/asset_loci.tsv"), show_col_types = FALSE)
regions <- loci %>%
  transmute(chr, start = pmax(1, pos - w), end = pos + w) %>%
  arrange(chr, start)
write_tsv(regions, file.path(nc$dir, "regions.txt"), col_names = FALSE)

# HapMap3 position map for 03_stream_ldsc.sh (Pan-UKB files have no rsIDs) ----
hm3 <- fread(cfg$ldsc$hm3)
pos <- rbindlist(lapply(1:22, function(k) {
  fread(paste0(cfg$ldsc$w_ld, k, ".l2.ldscore.gz"), select = c("CHR", "SNP", "BP"))
}))
map <- merge(pos, hm3, by = "SNP")[, .(chr = CHR, pos = BP, a1 = A1, a2 = A2, rsid = SNP)]
map <- map[!duplicated(paste(chr, pos)) & !duplicated(rsid)][order(chr, pos)]
fwrite(map, file.path(nc$dir, "hm3_map.tsv"), sep = "\t")
message(nrow(map), " HapMap3 SNPs in the position map")
