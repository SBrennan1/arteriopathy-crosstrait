#!/usr/bin/env bash
# Extracts the ASSET regions for each negative control from Pan-UKB with remote tabix
# (EUR columns; low-confidence variants removed).
# Usage: bash 10_negative_controls/02_extract_regions.sh [parallel jobs, default 8]
set -euo pipefail

JOBS=${1:-8}
DIR=${NC_DIR:-data/negative_controls}
DRAWS=$(pwd)/${NC_DRAWS:-config/negative_control_draws.tsv}
REF=$(pwd)/${NC_REF:-config/negative_control_reference.tsv}
inputs=("$DRAWS"); [[ -s $REF ]] && inputs+=("$REF")
REGIONS=$(pwd)/$DIR/regions.txt
export REGIONS

command -v tabix >/dev/null || { echo "tabix not found (htslib)" >&2; exit 1; }
[[ -s $DRAWS && -s $REGIONS ]] || { echo "run 01_select_controls.R first" >&2; exit 1; }
cd "$DIR"
mkdir -p regional tbi_cache logs

# distinct controls: id, url, index url
list=$(awk -F'\t' '
  FNR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
  !seen[$(c["id"])]++ { print $(c["id"]) "\t" $(c["aws_link"]) "\t" $(c["aws_link_tabix"]) }' "${inputs[@]}")

# The .tbi files sit in a separate folder from the data, so each index is fetched once
# into tbi_cache and passed to tabix explicitly (file##idx##index)
get_index() {
  local id=$1 tbi=$2 idx=$PWD/tbi_cache/$1.tsv.bgz.tbi
  [[ -s $idx ]] || curl -sfL -o "$idx" "$tbi" || { rm -f "$idx"; return 1; }
  echo "$idx"
}
export -f get_index

# region chromosome names have to match the file contigs ("1" vs "chr1")
read -r first_id first_url first_tbi <<< "$(printf '%s\n' "$list" | head -1)"
first_idx=$(get_index "$first_id" "$first_tbi") || { echo "could not fetch $first_tbi" >&2; exit 1; }
contigs=" $(tabix -l "$first_url##idx##$first_idx" | tr '\n' ' ')"
for ch in $(cut -f1 "$REGIONS" | sort -u); do
  [[ "$contigs" == *" $ch "* ]] || { echo "chromosome '$ch' not among file contigs:$contigs" >&2; exit 1; }
done
echo "$(wc -l < "$REGIONS") regions, $(printf '%s\n' "$list" | grep -c .) controls"

extract_one() {
  set -o pipefail
  local id=$1 url=$2 tbi=$3 out=regional/$1.tsv idx
  [[ -s $out ]] && return 0
  idx=$(get_index "$id" "$tbi") || { echo "FAIL $id: index not fetched" >&2; return 1; }

  # Header from the first compressed block. tabix -h does not return it when the index
  # was built with a skipped header line rather than a '#' one.
  local hdr
  hdr=$(curl -sfL -r 0-262143 "$url" | gzip -dc 2>/dev/null | head -1 || true)

  if ! tabix "$url##idx##$idx" -R "$REGIONS" | awk -F'\t' -v OFS='\t' -v H="$hdr" '
      BEGIN {
        n = split(H, h, "\t"); for (i = 1; i <= n; i++) c[h[i]] = i
        k = split("chr pos ref alt af_controls_EUR beta_EUR se_EUR neglog10_pval_EUR low_confidence_EUR", need, " ")
        for (j = 1; j <= k; j++) if (!(need[j] in c)) { print "missing column " need[j] > "/dev/stderr"; exit 2 }
        print "chr", "pos", "ref", "alt", "af_controls", "beta", "se", "nlp"
      }
      $1 == "chr" { next }
      tolower($(c["low_confidence_EUR"])) == "true" { next }
      $(c["se_EUR"]) == "NA" || $(c["se_EUR"]) == "" || $(c["se_EUR"]) + 0 == 0 { next }
      # overlapping regions return the same variant twice
      !seen[$(c["chr"]) ":" $(c["pos"]) ":" $(c["ref"]) ":" $(c["alt"])]++ {
        print $(c["chr"]), $(c["pos"]), $(c["ref"]), $(c["alt"]), $(c["af_controls_EUR"]),
              $(c["beta_EUR"]), $(c["se_EUR"]), $(c["neglog10_pval_EUR"])
      }' > "$out.tmp" 2> "logs/$id.regions.err"; then
    echo "FAIL $id (logs/$id.regions.err)" >&2; rm -f "$out.tmp"; return 1
  fi

  local n=$(( $(wc -l < "$out.tmp") - 1 ))
  [[ $n -gt 0 ]] || { echo "FAIL $id: no variants in the regions" >&2; rm -f "$out.tmp"; return 1; }
  mv "$out.tmp" "$out"
  echo "done $id: $n variants"
}
export -f extract_one

printf '%s\n' "$list" | xargs -P "$JOBS" -n 3 bash -c 'extract_one "$0" "$1" "$2"' \
  || echo "Some controls failed; rerun to retry them." >&2
