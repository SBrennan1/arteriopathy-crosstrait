# Shared genetic architecture across arteriopathies

Code for the cross-trait genetic analysis of five arteriopathies: cervical artery
dissection (CeAD), intracranial aneurysm (IA), spontaneous coronary artery dissection
(SCAD), aortic aneurysm and dissection (AAD; `AA` in the code) and fibromuscular
dysplasia (FMD).

> Brennan SO, et al. Shared polygenic architecture across arteriopathies: an integrative
> cross-trait analysis. *Circulation: Genomic and Precision Medicine* (2026).

## Analyses

| Folder | Analysis |
|---|---|
| `00_data` | Formatting of the GWAS summary statistics |
| `01_ldsc` | SNP heritability and genetic correlation (LDSC) |
| `02_multitrait` | ASSET, CPASSOC and pairwise MTAG |
| `03_coloc` | Cross-trait colocalization at the ASSET loci (coloc.abf, coloc.susie) |
| `04_gene_mapping` | Variant-to-gene mapping |
| `05_mr` | Proteome-wide and transcriptome-wide MR, pQTL replication |
| `06_qtl_coloc` | Colocalization of the MR associations |
| `07_downstream` | Genotype-Phenotype Map, blood pressure MR, enrichment, Open Targets, Table 1, drug target tiers |
| `08_figures` | Figures |
| `09_tables` | Supplementary tables |
| `10_negative_controls` | Pan-UK Biobank negative controls |
| `11_vascular_expression` | Cell-type specificity in single-cell atlases, AAA expression |
| `R` | Shared functions |

## Requirements

**R (>= 4.3):** ASSET, coloc, susieR, TwoSampleMR, ieugwasr, genetics.binaRies, qvalue,
clusterProfiler, ReactomePA, org.Hs.eg.db, gpmapr (mrcieu.r-universe.dev), rtracklayer, GenomicRanges, arrow,
reticulate, igraph, limma, GEOquery, data.table, vroom, tidyverse, yaml, httr, jsonlite,
openxlsx, ggplot2, patchwork, cowplot, ggrepel, ggrastr, eulerr, ComplexHeatmap,
locuszoomr, EnsDb.Hsapiens.v75.

**Python (>= 3.9):** numpy, scipy, pandas, h5py, anndata, pyyaml.

**Other tools:**
- [LDSC](https://github.com/bulik/ldsc) (`munge_sumstats.py`, `ldsc.py` on the `PATH`, or set in `config.yml`)
- [MTAG](https://github.com/JonJala/mtag) v1.0.8 (commit 9e17f3c), patched for Python 3:
  ```
  git clone https://github.com/JonJala/mtag external/mtag
  cd external/mtag && git checkout 9e17f3c && patch -p1 < ../../patches/mtag_python3.patch
  ```
- CPASSOC `FunctionSet.R` from the authors (Zhu et al. 2015), saved to `external/CPASSOC/`
- htslib (`tabix`) and `curl`, for the negative-control downloads
- [FUMA](https://fuma.ctglab.nl) v1.5.2 (SNP2GENE and GENE2FUNC, web)
- [Genotype-Phenotype Map](https://gpmap.opengwas.io) (GWAS upload, web)

API tokens are read from environment variables (for example in `~/.Renviron`):
`OPENGWAS_JWT` (OpenGWAS) and `LDLINK_TOKEN` (LDlink).

## Data

Input data are not included. See [`data/README.md`](data/README.md) for sources. Paths
are set in `config.yml`; machine-specific settings and the GPMap upload IDs go in
`config.local.yml` (see `config.local.yml.example`).

## Running

Run all scripts from the repository root. `run_all.R` lists the R scripts in order; the
steps outside R are:

1. **MTAG** after `00_data/02_munge_full.R`: `bash 02_multitrait/run_mtag.sh`
2. **FUMA SNP2GENE** on the MTAG statistics (`results/mtag/fuma_input/`) and on the
   genome-wide ASSET statistics (`results/asset/asset_all.tsv.gz`), and **GENE2FUNC** on
   the ASSET genes, with GRCh37, 1000 Genomes phase 3 EUR, lead P < 5e-8, r2 0.6 and 0.1,
   loci merged within 250 kb, MHC excluded and MAGMA with a 10 kb window. Save the
   outputs to the folders set under `mtag` and `enrichment` in `config.yml`.
3. **GPMap**: upload each formatted GWAS and set the upload IDs under `gpmap: guids`.
4. **Negative controls** after `10_negative_controls/01_select_controls.R`:
   `bash 10_negative_controls/02_extract_regions.sh` and
   `bash 10_negative_controls/03_stream_ldsc.sh`
5. **Single-cell summaries** after `11_vascular_expression/01_locus_groups.R`:
   `python 11_vascular_expression/02_expression_summary.py hca_artery` and
   `python 11_vascular_expression/02_expression_summary.py brain_aneurysm`

Summary counts from each stage are written to `results/counts.tsv`.

## Outputs

| Output | Script |
|---|---|
| Figure 2 | `08_figures/fig_rg_heatmap.R` |
| Figure 3 | `07_downstream/01_gpmap_bp.R` |
| Figure 4 | `08_figures/fig_mr_manhattan.R` |
| Figure 5 | `11_vascular_expression/05_figure.R` |
| Figure S1 | `08_figures/fig_miami.R` |
| Figure S2 | `08_figures/fig_locuszoom_novel.R` |
| Figure S3 | `10_negative_controls/06_figure.R` |
| Figure S4 | `08_figures/fig_gene2func.R` |
| Figure S5 | `08_figures/fig_method_euler.R` |
| Table 1 | `07_downstream/05_convergence_table.R` (`results/convergence/table1.tsv`) |
| Table 2 | `07_downstream/08_target_prioritization.R` (`results/prioritization/target_tiers_main.tsv`) |
| Supplementary Tables 1-31 | `09_tables/supplementary_tables.R` (`results/tables/supplementary_tables.xlsx`) |

Analyses use GRCh37; output tables also give GRCh38 positions.

## License

MIT. See [`LICENSE`](LICENSE).
