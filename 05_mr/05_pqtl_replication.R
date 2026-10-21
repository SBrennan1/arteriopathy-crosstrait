# Replication of FDR-significant pQTL MR associations on the same assay: deCODE
# (SomaScan) in Fenland and UKB-PPP (Olink) in SCALLOP. Cis-pQTLs (+/- 100 kb,
# P < 5e-8) clumped per protein (r2 < 0.001, 10 Mb). Replicated: same direction and
# P < 0.05 / number tested.

source("R/setup.R")
out <- res_dir("mr")
rp <- cfg$mr$replication

disc <- read_stage("mr", "mr_pqtl.tsv") %>%
  filter(qval < cfg$mr$fdr, platform %in% names(rp$cohort)) %>%
  transmute(platform, study = unlist(rp$cohort[platform]), protein = exposure, outcome,
            disc_b = b, disc_se = se, disc_q = qval)
message(nrow(disc), " discovery associations (", n_distinct(disc$protein), " proteins)")

in_window <- function(d, win) {
  d %>%
    inner_join(win %>% rename(w_chr = chr), by = "protein") %>%
    filter(chr == w_chr, pos >= start, pos <= end) %>%
    select(-w_chr, -start, -end)
}

# 1000 Genomes EUR variants, used to assign rsIDs to SCALLOP variants
bim <- fread(paste0(cfg$ld_ref, ".bim"), select = c(1, 2, 4, 5, 6),
             col.names = c("chr", "rsid", "pos", "a1", "a2"))
bim[, chr := suppressWarnings(as.integer(chr))]

# Fenland ----
fen_raw <- vroom(rp$fenland, show_col_types = FALSE) %>%
  mutate(across(where(is.character), ~ sub("^'", "", .x)))
fen <- fen_raw %>%
  filter(cis.trans == "cis", HGNC.symbol.protein %in% disc$protein[disc$study == "Fenland"]) %>%
  transmute(protein = HGNC.symbol.protein, rsid = rsID, chr = suppressWarnings(as.integer(Chr)),
            pos = as.integer(Position), effect_allele = toupper(EA), other_allele = toupper(NEA),
            eaf = as.numeric(EAF), beta = as.numeric(Effect), se = as.numeric(SE),
            pval = suppressWarnings(as.numeric(P.value)))
if (any(is.na(fen$pval))) {
  message(sum(is.na(fen$pval)), " Fenland P values below numeric range set to 1e-300")
  fen$pval[is.na(fen$pval)] <- 1e-300
}
fen <- fen %>%
  in_window(protein_regions(unique(fen$protein), kb = rp$window_kb)) %>%
  filter(str_detect(rsid, "^rs"), pval < rp$p) %>%
  arrange(pval) %>%
  distinct(protein, rsid, .keep_all = TRUE)

# SCALLOP ----
if (!dir.exists(rp$scallop$dir)) stop("SCALLOP directory not found: ", rp$scallop$dir, call. = FALSE)
scallop_file <- function(g) file.path(rp$scallop$dir, paste0(rp$scallop$file_names[[g]] %||% g, ".txt.gz"))
olink_prot <- unique(disc$protein[disc$study == "SCALLOP"])
on_scallop <- olink_prot[file.exists(map_chr(olink_prot, scallop_file))]
message("Olink discoveries with a SCALLOP file: ", paste(on_scallop, collapse = ", "))
message("Not on the SCALLOP panel: ", paste(setdiff(olink_prot, on_scallop), collapse = ", "))

win_o <- protein_regions(on_scallop, kb = rp$window_kb)
sc <- map_dfr(on_scallop, function(g) {
  d <- vroom(scallop_file(g), show_col_types = FALSE)
  if (!"Freq1" %in% names(d)) d$Freq1 <- NA_real_
  m <- str_split_fixed(d$MarkerName, ":", 3)
  tibble(protein = g, chr = suppressWarnings(as.integer(m[, 1])), pos = as.integer(m[, 2]),
         effect_allele = toupper(d$Allele1), other_allele = toupper(d$Allele2),
         eaf = as.numeric(d$Freq1), beta = d$Effect, se = d$StdErr, pval = d$`P-value`) %>%
    filter(pval < rp$p) %>%
    in_window(filter(win_o, protein == g))
}) %>%
  inner_join(as_tibble(bim), by = c("chr", "pos"), relationship = "many-to-many") %>%
  filter((effect_allele == a1 & other_allele == a2) | (effect_allele == a2 & other_allele == a1)) %>%
  select(-a1, -a2) %>%
  distinct(protein, rsid, .keep_all = TRUE)

# Clump each protein separately ----
clump_each <- function(d) {
  map_dfr(split(d, d$protein), function(x) {
    if (nrow(x) <= 1) return(x)
    keep <- clump_local(x, p = 1, r2 = rp$clump_r2, kb = rp$clump_kb)
    filter(x, rsid %in% keep$rsid)
  })
}
ivs <- bind_rows(Fenland = clump_each(fen), SCALLOP = clump_each(sc), .id = "study") %>%
  rename(exposure = protein)
print(count(ivs, study, exposure, name = "instruments"), n = Inf)

# MR ----
outcomes <- map(set_names(unique(disc$outcome)), read_mr_outcome, rsids = unique(ivs$rsid))
rep <- disc %>%
  distinct(study, outcome) %>%
  pmap_dfr(function(study, outcome) {
    pr <- disc$protein[disc$study == study & disc$outcome == outcome]
    x <- filter(ivs, study == !!study, exposure %in% pr)
    if (nrow(x) == 0) return(NULL)
    r <- run_mr(x, outcomes[[outcome]], outcome)
    if (is.null(r)) NULL else mutate(r, study = study)
  })

res <- disc %>%
  left_join(rep %>% transmute(study, protein = exposure, outcome, rep_method = method, rep_nsnp = nsnp,
                              rep_b = b, rep_se = se, rep_p = pval, rep_snps = snps),
            by = c("study", "protein", "outcome")) %>%
  mutate(tested = !is.na(rep_b),
         not_tested = case_when(
           tested ~ NA_character_,
           study == "SCALLOP" & !protein %in% on_scallop ~ "Not measured in SCALLOP",
           !paste(study, protein) %in% paste(ivs$study, ivs$exposure) ~ "No cis-pQTL instrument",
           TRUE ~ "Instruments not in the outcome GWAS"))
threshold <- 0.05 / sum(res$tested)
res <- res %>%
  mutate(concordant = if_else(tested, sign(rep_b) == sign(disc_b), NA),
         replicated = if_else(tested, concordant & rep_p < threshold, NA),
         replicated_nominal = if_else(tested, concordant & rep_p < 0.05, NA)) %>%
  arrange(desc(tested), rep_p, platform, outcome, protein)
write_tsv(res, file.path(out, "pqtl_replication.tsv"))

t <- filter(res, tested)
message(sprintf("%d of %d discovery associations tested; %d same direction, %d replicated at P < %.2g, %d at P < 0.05",
                nrow(t), nrow(res), sum(t$concordant), sum(t$replicated), threshold, sum(t$replicated_nominal)))
print(t %>% transmute(study, protein, outcome, rep_nsnp, disc_b = round(disc_b, 3), rep_b = round(rep_b, 3),
                      rep_p = signif(rep_p, 3), concordant, replicated), n = Inf)
print(count(res, study, not_tested))

write_counts("pqtl_replication", c(
  discovery = nrow(res),
  tested = nrow(t),
  tested_fenland = sum(t$study == "Fenland"),
  tested_scallop = sum(t$study == "SCALLOP"),
  concordant = sum(t$concordant),
  replicated_bonferroni = sum(t$replicated),
  replicated_nominal = sum(t$replicated_nominal),
  bonferroni_threshold = signif(threshold, 3)
))
