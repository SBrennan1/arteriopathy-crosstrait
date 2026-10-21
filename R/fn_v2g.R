# Variant-to-gene mapping shared by the ASSET, CPASSOC and MTAG analyses.
# Tiers, in order: coding variant (Ensembl VEP, GRCh37), eQTL colocalisation
# with GTEx v10 arterial tissues (coloc.abf, PP.H4 > 0.8), Open Targets E2G
# interval score > 0.5, nearest protein-coding gene. Curated overrides from
# config/v2g_overrides.tsv are applied afterwards and kept in separate columns.

coding_consequences <- c(
  "missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
  "splice_acceptor_variant", "splice_donor_variant", "inframe_insertion", "inframe_deletion"
)

# Most severe consequence per rsID, with the gene of the first protein-coding
# transcript consequence
vep_annotate <- function(rsids, url = cfg$v2g$vep_url, batch_size = 200) {
  rsids <- unique(na.omit(rsids))
  empty <- tibble(rsid = character(), vep_consequence = character(), vep_gene = character(),
                  vep_gene_id = character(), is_coding = logical())
  if (length(rsids) == 0) return(empty)
  raw <- map(split(rsids, ceiling(seq_along(rsids) / batch_size)), function(ids) {
    r <- httr::POST(url, httr::content_type_json(), httr::accept_json(),
                    body = jsonlite::toJSON(list(ids = ids)), encode = "raw")
    Sys.sleep(0.5)
    if (httr::status_code(r) != 200) {
      warning("VEP request failed with status ", httr::status_code(r))
      return(NULL)
    }
    jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), flatten = TRUE)
  })
  vep <- bind_rows(raw)
  if (nrow(vep) == 0) return(empty)
  vep <- distinct(vep, id, .keep_all = TRUE)

  tcs <- vep$transcript_consequences %||% vector("list", nrow(vep))
  first_gene <- map_dfr(tcs, function(tc) {
    if (!is.data.frame(tc) || nrow(tc) == 0) {
      return(tibble(vep_gene = NA_character_, vep_gene_id = NA_character_))
    }
    pc <- which(tc$biotype == "protein_coding")
    if (length(pc)) tc <- tc[pc, , drop = FALSE]
    tibble(vep_gene = tc$gene_symbol[1] %||% NA_character_,
           vep_gene_id = tc$gene_id[1] %||% NA_character_)
  })
  tibble(rsid = vep$id, vep_consequence = vep$most_severe_consequence) %>%
    bind_cols(first_gene) %>%
    mutate(is_coding = vep_consequence %in% coding_consequences)
}

# GTEx gene IDs to symbols via the protein-coding GTF genes; anything else
# (lncRNAs etc.) keeps its versioned Ensembl ID
ensg_symbol <- function(ids, genes) {
  g <- distinct(genes, gene_id, .keep_all = TRUE)
  coalesce(g$gene_name[match(sub("\\..*", "", ids), g$gene_id)], ids)
}

query_e2g <- function(variant_id, url = cfg$v2g$opentargets_url) {
  q <- sprintf('{
    variant(variantId: "%s") {
      id rsIds
      intervals {
        count
        rows {
          target { id approvedSymbol }
          score datasourceId biosampleName intervalType distanceToTss
        }
      }
    }
  }', variant_id)
  r <- httr::POST(url, httr::content_type_json(),
                  body = jsonlite::toJSON(list(query = q), auto_unbox = TRUE), encode = "raw")
  if (httr::status_code(r) != 200) return(NULL)
  jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), flatten = TRUE)$data$variant
}

# Highest-scoring Open Targets E2G interval per variant. d: rsid, chr_hg38,
# pos_hg38, a1, a2. Variant IDs are chr_pos_ref_alt on GRCh38 and the reference
# allele is not known here, so both allele orders are tried.
e2g_annotate <- function(d) {
  d <- d %>% filter(!is.na(pos_hg38)) %>% distinct(rsid, .keep_all = TRUE)
  rows <- pmap_dfr(d %>% select(rsid, chr_hg38, pos_hg38, a1, a2), function(rsid, chr_hg38, pos_hg38, a1, a2) {
    res <- NULL
    for (id in sprintf("%s_%d_%s_%s", chr_hg38, as.integer(pos_hg38), c(a2, a1), c(a1, a2))) {
      res <- tryCatch(query_e2g(id), error = function(e) NULL)
      if (!is.null(res)) break
    }
    Sys.sleep(0.3)
    if (!isTRUE(res$intervals$count > 0)) return(NULL)
    as_tibble(res$intervals$rows) %>% mutate(rsid = rsid)
  })
  if (nrow(rows) == 0) {
    return(tibble(rsid = character(), e2g_gene = character(), e2g_score = numeric(), e2g_tissue = character()))
  }
  rows %>%
    group_by(rsid) %>%
    summarise(e2g_gene = target.approvedSymbol[which.max(score)], e2g_score = max(score),
              e2g_tissue = biosampleName[which.max(score)], .groups = "drop")
}

# GTEx eQTLs within w of any centre (GRCh38), one parquet per tissue x chromosome
load_gtex_windows <- function(centres, w, gtex = cfg$gtex_v10) {
  todo <- expand_grid(tissue = names(gtex$tissues), chr = unique(as.character(centres$chr)))
  files <- file.path(gtex$dir, sprintf("%s.v10.allpairs.chr%s.parquet", todo$tissue, todo$chr))
  if (any(!file.exists(files))) {
    warning("GTEx files not found, skipped: ", paste(basename(files[!file.exists(files)]), collapse = ", "))
  }
  map_dfr(which(file.exists(files)), function(i) {
    ctr <- centres$pos_hg38[as.character(centres$chr) == todo$chr[i]]
    d <- arrow::read_parquet(files[i], col_select = c("gene_id", "variant_id", "slope", "slope_se")) %>%
      mutate(eqtl_pos = as.integer(str_extract(variant_id, "(?<=_)\\d+(?=_)")))
    keep <- Reduce(`|`, lapply(ctr, function(x) d$eqtl_pos >= x - w & d$eqtl_pos <= x + w))
    d[which(keep), ] %>%
      transmute(tissue = todo$tissue[i], chr = todo$chr[i], gene_id, eqtl_pos, slope, slope_se)
  })
}

coloc_abf_quiet <- function(...) {
  fit <- NULL
  invisible(capture.output(fit <- tryCatch(suppressMessages(coloc::coloc.abf(...)), error = function(e) NULL)))
  fit
}

# coloc.abf between each GWAS trait and every GTEx gene in the window around each
# locus centre. centres: locus_rsid, chr, pos (GRCh37 centre), traits (comma-
# separated GWAS traits). GWAS and eQTL variants are matched on GRCh38 position.
gtex_coloc <- function(centres, p = cfg$v2g, gtex = cfg$gtex_v10) {
  w <- p$window_kb * 1e3
  centres <- centres %>%
    select(locus_rsid, chr, pos, traits) %>%
    add_hg38() %>%
    filter(!is.na(pos_hg38))
  jobs <- centres %>%
    separate_rows(traits, sep = ",\\s*") %>%
    rename(trait = traits) %>%
    filter(trait %in% names(cfg$traits))

  gwas <- map_dfr(unique(jobs$trait), function(t) {
    by_chr <- read_univar(t, cols = c("rsid", "chr", "pos", "beta", "se", "N")) %>%
      filter(!is.na(beta), !is.na(se), se > 0) %>%
      split(.$chr)
    j <- filter(jobs, trait == t)
    map_dfr(seq_len(nrow(j)), function(i) {
      d <- by_chr[[as.character(j$chr[i])]]
      if (is.null(d)) return(NULL)
      d %>%
        filter(pos >= j$pos[i] - w, pos <= j$pos[i] + w) %>%
        distinct(rsid, .keep_all = TRUE) %>%
        mutate(locus_rsid = j$locus_rsid[i], trait = t)
    })
  }) %>%
    add_hg38() %>%
    filter(!is.na(pos_hg38), chr_hg38 == as.character(chr))

  eqtl <- load_gtex_windows(centres, w, gtex) %>%
    filter(!is.na(slope), !is.na(slope_se), slope_se > 0)

  res <- pmap_dfr(jobs %>% select(locus_rsid, trait, chr, pos_hg38), function(locus_rsid, trait, chr, pos_hg38) {
    g <- gwas[gwas$locus_rsid == locus_rsid & gwas$trait == trait, ]
    if (nrow(g) < p$min_gwas_snps) return(NULL)
    n_gwas <- suppressWarnings(max(g$N, na.rm = TRUE))
    if (!is.finite(n_gwas)) n_gwas <- cfg$traits[[trait]]$ncase + cfg$traits[[trait]]$ncontrol
    e <- eqtl[eqtl$chr == as.character(chr) & eqtl$eqtl_pos >= pos_hg38 - w & eqtl$eqtl_pos <= pos_hg38 + w, ]
    if (nrow(e) == 0) return(NULL)

    map_dfr(split(e, list(e$tissue, e$gene_id), drop = TRUE), function(eg) {
      eg <- distinct(eg, eqtl_pos, .keep_all = TRUE)
      if (nrow(eg) < p$min_eqtl_snps) return(NULL)
      m <- inner_join(g %>% select(rsid, pos_hg38, beta, se),
                      eg %>% select(eqtl_pos, slope, slope_se),
                      by = c(pos_hg38 = "eqtl_pos")) %>%
        distinct(pos_hg38, .keep_all = TRUE)
      if (nrow(m) < p$min_eqtl_snps) return(NULL)
      fit <- coloc_abf_quiet(
        dataset1 = list(beta = m$beta, varbeta = m$se^2, snp = m$rsid, type = "cc", s = 0.5, N = n_gwas),
        dataset2 = list(beta = m$slope, varbeta = m$slope_se^2, snp = m$rsid, type = "quant",
                        N = gtex$tissues[[eg$tissue[1]]], sdY = 1),
        p1 = p$p1, p2 = p$p2, p12 = p$p12)
      if (is.null(fit)) return(NULL)
      tibble(locus_rsid = locus_rsid, gwas_trait = trait, tissue = eg$tissue[1],
             gene_id = eg$gene_id[1], n_snps = nrow(m),
             as_tibble(as.list(fit$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])))
    })
  })
  if (nrow(res) == 0) {
    res <- tibble(locus_rsid = character(), gwas_trait = character(), tissue = character(),
                  gene_id = character(), n_snps = integer(), PP.H0.abf = numeric(),
                  PP.H1.abf = numeric(), PP.H2.abf = numeric(), PP.H3.abf = numeric(),
                  PP.H4.abf = numeric())
  }
  res
}

# Best eQTL colocalisation per locus above the PP.H4 threshold
best_eqtl_coloc <- function(coloc_all, genes, min_h4 = cfg$v2g$pp_h4) {
  coloc_all %>%
    filter(PP.H4.abf > min_h4) %>%
    group_by(locus_rsid) %>%
    slice_max(PP.H4.abf, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(locus_rsid, coloc_gene = ensg_symbol(gene_id, genes), coloc_gene_id = gene_id,
              coloc_tissue = tissue, coloc_pp_h4 = PP.H4.abf, coloc_trait = gwas_trait,
              coloc_n_snps = n_snps)
}

# Needs coding_gene, coloc_gene, e2g_gene, e2g_score and nearest_gene
assign_tiers <- function(d, min_e2g = cfg$v2g$e2g_min_score) {
  d %>%
    mutate(
      .e2g = !is.na(e2g_gene) & coalesce(e2g_score > min_e2g, FALSE),
      gene_tiered = case_when(
        !is.na(coding_gene) ~ coding_gene,
        !is.na(coloc_gene) ~ coloc_gene,
        .e2g ~ e2g_gene,
        TRUE ~ nearest_gene
      ),
      evidence_tiered = case_when(
        !is.na(coding_gene) ~ "Coding variant",
        !is.na(coloc_gene) ~ "eQTL colocalization",
        .e2g ~ "Chromatin interaction",
        TRUE ~ "Nearest gene"
      )
    ) %>%
    select(-.e2g)
}

# Curated assignments keyed by lead rsID. An empty gene or evidence in the file
# keeps the tiered value (the row then only adds a note).
apply_v2g_overrides <- function(d, analysis, key = "rsid", path = cfg$v2g$overrides) {
  ov <- read_tsv(path, col_types = cols(.default = "c"), na = "") %>%
    filter(analysis == !!analysis) %>%
    transmute(.key = rsid, .gene = gene, .evidence = evidence, override_note = note)
  missing <- setdiff(ov$.key, d[[key]])
  if (length(missing)) {
    warning(analysis, " overrides not applied, rsID not among the leads: ", paste(missing, collapse = ", "))
  }
  d %>%
    left_join(ov, by = setNames(".key", key)) %>%
    mutate(gene = coalesce(.gene, gene_tiered),
           evidence = coalesce(.evidence, evidence_tiered),
           overridden = !is.na(.gene) & !coalesce(.gene == gene_tiered, FALSE)) %>%
    select(-.gene, -.evidence) %>%
    relocate(override_note, .after = overridden)
}

# Full tiered mapping for single lead variants (CPASSOC, MTAG).
# leads: rsid, chr, pos, a1, a2, traits (comma-separated GWAS traits for the eQTL coloc)
v2g_leads <- function(leads, analysis) {
  genes <- load_genes()
  leads <- leads %>% add_hg38()
  vep <- vep_annotate(leads$rsid)
  coloc_all <- gtex_coloc(leads %>% transmute(locus_rsid = rsid, chr, pos, traits))
  e2g <- e2g_annotate(leads)

  v2g <- leads %>%
    bind_cols(nearest_gene(leads$chr, leads$pos, genes)) %>%
    left_join(vep, by = "rsid") %>%
    mutate(coding_gene = if_else(coalesce(is_coding, FALSE), vep_gene, NA_character_)) %>%
    left_join(best_eqtl_coloc(coloc_all, genes), by = c(rsid = "locus_rsid")) %>%
    left_join(e2g, by = "rsid") %>%
    assign_tiers() %>%
    apply_v2g_overrides(analysis) %>%
    select(rsid, chr, pos, chr_hg38, pos_hg38, traits, vep_consequence, coding_gene,
           starts_with("coloc_"), starts_with("e2g_"), nearest_gene, nearest_distance,
           gene_tiered, evidence_tiered, gene, evidence, overridden, override_note)

  list(v2g = v2g, coloc = coloc_all %>% mutate(gene = ensg_symbol(gene_id, genes), .after = gene_id))
}
