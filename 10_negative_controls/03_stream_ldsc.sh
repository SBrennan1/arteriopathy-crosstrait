#!/usr/bin/env bash
# Streams each negative control used for genetic correlation from Pan-UKB, keeping
# HapMap3 SNPs (ldsc/<id>.tsv.gz) and variants with P < 5e-8 (gws/<id>.tsv).
# Usage: bash 10_negative_controls/03_stream_ldsc.sh [parallel jobs, default 4]
set -euo pipefail

JOBS=${1:-4}
DIR=${NC_DIR:-data/negative_controls}
DRAWS=$(pwd)/${NC_DRAWS:-config/negative_control_draws.tsv}
REF=$(pwd)/${NC_REF:-config/negative_control_reference.tsv}
inputs=("$DRAWS"); [[ -s $REF ]] && inputs+=("$REF")

[[ -s $DRAWS && -s $DIR/hm3_map.tsv ]] || { echo "run 01_select_controls.R first" >&2; exit 1; }
cd "$DIR"
mkdir -p ldsc gws logs

# distinct LDSC-subset controls: id, url, total N
list=$(awk -F'\t' '
  FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
  $(c["ldsc_subset"]) == "TRUE" && !seen[$(c["id"])]++ {
    print $(c["id"]) "\t" $(c["aws_link"]) "\t" ($(c["n_cases_EUR"]) + $(c["n_controls_EUR"]))
  }' "${inputs[@]}")
n_traits=$(printf '%s\n' "$list" | grep -c . || true)
[[ $n_traits -gt 0 ]] || { echo "no LDSC-subset controls in $DRAWS" >&2; exit 1; }

first_url=$(printf '%s\n' "$list" | head -1 | cut -f2)
curl -sIL "$first_url" | awk -v n="$n_traits" 'tolower($1) == "content-length:" {
  gb = $2 / 1e9; printf "First file %.1f GB, roughly %.0f GB for %d controls\n", gb, gb * n, n }'

stream_one() {
  set -o pipefail
  local id=$1 url=$2 n=$3
  [[ -s ldsc/$id.tsv.gz && -s gws/$id.tsv ]] && return 0

  if ! curl -sfL "$url" | gzip -dc | awk -F'\t' -v OFS='\t' -v N="$n" \
        -v H="ldsc/$id.tmp.gz" -v G="gws/$id.tsv.tmp" '
      BEGIN { comp["A"] = "T"; comp["T"] = "A"; comp["C"] = "G"; comp["G"] = "C"; zh = "gzip -c > " H }
      # hm3_map.tsv: chr pos a1 a2 rsid
      NR == FNR { if (FNR > 1) { rs[$1 ":" $2] = $5; al[$1 ":" $2] = $3 $4 }; next }
      FNR == 1 {
        for (i = 1; i <= NF; i++) c[$i] = i
        k = split("chr pos ref alt af_controls_EUR beta_EUR se_EUR neglog10_pval_EUR low_confidence_EUR", need, " ")
        for (j = 1; j <= k; j++) if (!(need[j] in c)) { print "missing column " need[j] > "/dev/stderr"; exit 2 }
        print "rsid", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "N" | zh
        print "chr", "pos", "ref", "alt", "af_controls", "beta", "se", "nlp" > G
        next
      }
      tolower($(c["low_confidence_EUR"])) == "true" { next }
      {
        b = $(c["beta_EUR"]); s = $(c["se_EUR"]); nlp = $(c["neglog10_pval_EUR"])
        if (b == "NA" || s == "NA" || b == "" || s == "" || s + 0 == 0) next
        ref = $(c["ref"]); alt = $(c["alt"]); key = $(c["chr"]) ":" $(c["pos"])
        # keep the HapMap3 allele pair only (multi-allelic sites), either strand
        if (key in rs && (al[key] == ref alt || al[key] == alt ref ||
                          al[key] == comp[ref] comp[alt] || al[key] == comp[alt] comp[ref])) {
          p = (nlp + 0 > 300) ? 1e-300 : exp(-nlp * log(10))
          print rs[key], alt, ref, $(c["af_controls_EUR"]), b, s, p, N | zh
        }
        if (nlp + 0 > 7.30103)
          print $(c["chr"]), $(c["pos"]), ref, alt, $(c["af_controls_EUR"]), b, s, nlp > G
      }
      END { close(zh) }' hm3_map.tsv - 2> "logs/$id.ldsc.err"; then
    echo "FAIL $id (logs/$id.ldsc.err)" >&2; return 1
  fi

  local nh=$(( $(gzip -dc "ldsc/$id.tmp.gz" | wc -l) - 1 ))
  local ng=$(( $(wc -l < "gws/$id.tsv.tmp") - 1 ))
  if [[ $nh -lt 500000 ]]; then
    echo "FAIL $id: only $nh HapMap3 rows, left as .tmp" >&2; return 1
  fi
  mv "ldsc/$id.tmp.gz" "ldsc/$id.tsv.gz"
  mv "gws/$id.tsv.tmp" "gws/$id.tsv"
  echo "done $id: $nh HapMap3 rows, $ng variants at P < 5e-8"
}
export -f stream_one

printf '%s\n' "$list" | xargs -P "$JOBS" -n 3 bash -c 'stream_one "$0" "$1" "$2"' \
  || echo "Some controls failed; rerun to retry them." >&2
