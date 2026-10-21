#!/usr/bin/env bash
# Pairwise MTAG (v1.0.8 with patches/mtag_python3.patch).
# Run from the repository root after 00_data/02_munge_full.R.
set -euo pipefail

MTAG=${MTAG:-external/mtag/mtag.py}
PYTHON=${PYTHON:-python}
IN=results/munged_full
OUT=data/mtag
LDREF=${LDREF:-data/ref/mtag_ld/1000G_Phase3_weights_hm3_no_MHC.}

mkdir -p "$OUT"

run_pair() {
  "$PYTHON" "$MTAG" \
    --sumstats "$IN/$1.sumstats.gz,$IN/$2.sumstats.gz" \
    --out "$OUT/$1_$2" \
    --snp_name SNP --z_name Z --n_name N --eaf_name FRQ \
    --a1_name A1 --a2_name A2 --bpos_name BP --no_chr_data \
    --ld_ref_panel "$LDREF" \
    --fdr --stream_stdout
}

run_pair SCAD CeAD
run_pair AA CeAD
run_pair AA IA
