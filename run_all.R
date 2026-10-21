# Runs the R scripts in order. Steps outside R are run separately where noted (see
# README): MTAG, FUMA, GPMap upload, the negative-control downloads and the single-cell
# expression summaries.

scripts <- c(
  "00_data/01_prepare_sumstats.R",
  "00_data/02_munge_full.R",
  "01_ldsc/01_ldsc_rg.R",
  "01_ldsc/02_h2_liability.R",
  "02_multitrait/01_asset.R",
  "02_multitrait/02_cpassoc.R",
  # bash 02_multitrait/run_mtag.sh
  "02_multitrait/03_mtag_fuma_prep.R",
  # FUMA SNP2GENE on the MTAG and ASSET statistics
  "02_multitrait/04_mtag_leads.R",
  "02_multitrait/05_method_overlap.R",
  "02_multitrait/06_mtag_parameters.R",
  "02_multitrait/07_intercept_sensitivity.R",
  "03_coloc/01_asset_coloc_abf.R",
  "03_coloc/02_asset_susie_coloc.R",
  "03_coloc/03_cross_trait_coloc.R",
  "04_gene_mapping/01_v2g_asset.R",
  "04_gene_mapping/02_v2g_mtag.R",
  "04_gene_mapping/03_v2g_cpassoc.R",
  "05_mr/00_eqtl_instruments.R",
  "05_mr/01_instruments_olink.R",
  "05_mr/02_mr_pqtl.R",
  "05_mr/03_mr_eqtl.R",
  "05_mr/04_fdr.R",
  "05_mr/05_pqtl_replication.R",
  "06_qtl_coloc/00_eqtl_subsets.R",
  "06_qtl_coloc/01_pqtl_coloc.R",
  "06_qtl_coloc/02_eqtl_coloc.R",
  "06_qtl_coloc/03_mr_coloc_tables.R",
  "07_downstream/01_gpmap_bp.R",
  "07_downstream/02_enrichment.R",
  "07_downstream/03_magma.R",
  "07_downstream/04_drug_targets.R",
  "07_downstream/05_convergence_table.R",
  "07_downstream/06_cytobands.R",
  "07_downstream/07_pqtl_protein_altering.R",
  "10_negative_controls/01_select_controls.R",
  # bash 10_negative_controls/02_extract_regions.sh; bash 10_negative_controls/03_stream_ldsc.sh
  "10_negative_controls/04_ldsc_controls.R",
  "10_negative_controls/05_regional_null.R",
  "11_vascular_expression/01_locus_groups.R",
  # python 11_vascular_expression/02_expression_summary.py hca_artery (and brain_aneurysm)
  "11_vascular_expression/03_celltype_specificity.R",
  "11_vascular_expression/04_aneurysm_datasets.R",
  "07_downstream/08_target_prioritization.R",
  "08_figures/fig_rg_heatmap.R",
  "08_figures/fig_mr_manhattan.R",
  "08_figures/fig_miami.R",
  "08_figures/fig_locuszoom_novel.R",
  "08_figures/fig_gene2func.R",
  "08_figures/fig_method_euler.R",
  "10_negative_controls/06_figure.R",
  "11_vascular_expression/05_figure.R",
  "09_tables/supplementary_tables.R"
)

for (s in scripts) {
  message(s)
  if (system2("Rscript", s) != 0) stop("failed: ", s)
}
