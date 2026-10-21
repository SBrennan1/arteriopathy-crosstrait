# Input data

Input data are not distributed with this repository. Download them from the sources
below and set the paths in `config.yml` (the defaults are shown). Files marked
*downloaded by the code* are fetched by the scripts.

## GWAS summary statistics (European ancestry)

| Trait | Study | Source | Path |
|---|---|---|---|
| SCAD | Adlam et al. 2023, *Nat Genet* | GWAS Catalog GCST90245878 (GRCh37) | `raw/GCST90245878_buildGRCh37.tsv.gz` |
| CeAD | Debette et al. 2015, *Nat Genet* (CADISP) | On request from the CADISP investigators | `raw/cadisp.CEADALL.Debette2015_1000G.EUR.GWAS.txt.gz` |
| FMD | Georges et al. 2021, *Nat Commun* | GWAS Catalog GCST90026612 (harmonised) | `raw/GCST90026612.h.tsv.gz` |
| IA | Bakker et al. 2020, *Nat Genet* | Stage 1 excluding UK Biobank, from the study authors | `raw/IA.GWAS.BakkerMK.2020.sumstats.Stage_1_excludingUKBB.txt.gz` |
| AAD | Million Veteran Program, phecode 442.1 | GWAS Catalog GCST90475993 | `raw/GCST90475993.tsv.gz` |
| CAD | Aragam et al. 2022, *Nat Genet* | GWAS Catalog GCST90132314 (GRCh37) | `ref/GCST90132314_buildGRCh37.tsv` |

Blood pressure GWAS (OpenGWAS `ieu-b-5138`, `ieu-b-5139`, `ieu-b-5140`) are
*downloaded by the code*.

## Reference data

| Data | Source | Path |
|---|---|---|
| LDSC LD scores, weights and HapMap3 SNP list (1000 Genomes phase 3 EUR, no MHC) | LDSC reference files | `ref/1000G_Phase3_ldscores/`, `ref/1000G_Phase3_weights_hm3_no_MHC/`, `ref/hm3_no_MHC_3col.txt` |
| 1000 Genomes phase 3 EUR PLINK files, genome-wide and per chromosome | http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz | `ref/1kg_eur/`, `ref/1kg_eur_chr/` |
| UK Biobank LD matrices (Weissbrod et al. 2020), chunks covering the ASSET loci | `s3://broad-alkesgroup-ukbb-ld/UKBB_LD/` | `ref/ukbb_ld/` |
| GTEx v10 eQTL all-pairs (parquet), aorta, coronary and tibial artery, whole blood, cultured fibroblasts | GTEx Portal (`gs://gtex-resources/GTEx_Analysis_v10_QTLs/`) | `ref/gtex_v10/` |
| GTEx v10 gene median TPM | GTEx Portal | `ref/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz` |
| GENCODE v39 annotation (GRCh38) | GENCODE | `ref/gencode.v39.annotation.gtf.gz` |
| Ensembl GRCh37 release 87 GTF | Ensembl | `ref/Homo_sapiens.GRCh37.87.gtf.gz` |
| liftOver chains hg19ToHg38, hg38ToHg19 | UCSC | `ref/` |
| GWAS Catalog associations | GWAS Catalog | `ref/gwas_catalog_associations.tsv` |

## QTL data

| Data | Source | Path |
|---|---|---|
| UKB-PPP (Olink) cis-pQTL table | Sun et al. 2023, *Nature*, Supplementary Table | `mr/Olink.csv` |
| UKB-PPP discovery summary statistics and rsID maps | Synapse (UKB-PPP) | `pqtl/ukb_ppp/`, `pqtl/ukb_ppp_rsid_map/` |
| deCODE (SomaScan) cis-pQTL instruments, GRCh37 | Ferkingstad et al. 2021, *Nat Genet* | `mr/pqtls_soma_37.csv` |
| deCODE regional summary statistics per protein, lifted to GRCh37 | https://www.decode.com/summarydata/ | `decode/37/{PROTEIN}.csv` |
| GTEx v10 lead cis-eQTL per gene (P < 5e-8), five tissues | GTEx Portal | `mr/eqtl/{tissue}_all_genes_QCd_mr_instruments.csv` |
| Fenland (SomaScan) pQTLs | Pietzner et al. 2021, *Science* | `mr/replication/Fenland.csv` |
| SCALLOP CVD-I (Olink) per-protein summary statistics | Folkersen et al. 2020, *Nat Metab* | `mr/replication/{protein}.txt.gz` |

## Expression data

| Data | Source | Path |
|---|---|---|
| Human Cell Atlas arterial atlas (Zhao et al. 2025, *Cell Genomics*) | CELLxGENE dataset 43245158-5ae1-4e71-a9a6-67eef49c26bc | `vascular/hca_artery.h5ad` |
| Brain aneurysm atlas (Wang et al. 2026, *Nat Neurosci*) | https://brain-aneur.cells.ucsc.edu (`counts_*`, `meta.tsv`) | `vascular/brain_aneurysm/` |
| Abdominal aortic aneurysm vs donor aorta | GEO GSE57691, *downloaded by the code* | `vascular/gse57691/` |

## Negative controls

Pan-UK Biobank summary statistics for the control phenotypes in
`config/negative_control_draws.tsv` and `config/negative_control_reference.tsv` are
*downloaded by the code* (`10_negative_controls/`) to `negative_controls/`.

## Outputs of external tools

| Data | Path |
|---|---|
| Pairwise MTAG output and logs (`02_multitrait/run_mtag.sh`) | `mtag/` |
| FUMA SNP2GENE output for each MTAG target trait | `fuma/{trait}/` |
| FUMA SNP2GENE (MAGMA) output for the ASSET statistics | `fuma/ASSET/` |
| FUMA GENE2FUNC output for the ASSET genes | `fuma/gene2func_ASSET/` |

## Tables prepared outside the pipeline

Supplementary Tables 11 and 12 (cross-ancestry evaluation in East Asian IA and African
ancestry AAD GWAS) are read from `cross_ancestry/`. Supplementary Table 1 is
`config/study_characteristics.tsv`.
