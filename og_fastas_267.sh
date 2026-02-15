#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Subset ML-selected orthogroup FASTAs (1154 species) down to your 267 species
# -----------------------------

OGDIR="data/raw/y1000p_orthofinder/Orthogroup_Sequences"
OUTDIR="data/processed/03_ml_orthogroups/05_og_fastas_267"
ALOG="data/processed/03_ml_orthogroups/logs"
LIST="data/processed/03_ml_orthogroups/01_og_lists/orthogroup_list.txt"
SP="data/processed/species_267.txt"

mkdir -p "$OUTDIR" "$ALOG"

# Logs
MISSING_LOG="${ALOG}/missing_og_fastas.log"
COUNTS_TSV="${ALOG}/og_counts_267.tsv"
ZERO_LOG="${ALOG}/zero_after_filter.log"

: > "$MISSING_LOG"
: > "$COUNTS_TSV"
: > "$ZERO_LOG"

# Safety checks
if [ ! -s "$LIST" ]; then
  echo "[ERROR] OG list not found or empty: $LIST" >&2
  exit 1
fi

if [ ! -s "$SP" ]; then
  echo "[ERROR] species_267.txt not found or empty: $SP" >&2
  exit 1
fi

if [ ! -d "$OGDIR" ]; then
  echo "[ERROR] OGDIR does not exist: $OGDIR" >&2
  exit 1
fi

# Optional: quick dependency check
command -v seqkit >/dev/null 2>&1 || { echo "[ERROR] seqkit not found in PATH" >&2; exit 1; }

# Header for counts
echo -e "OG\tNseq_267" >> "$COUNTS_TSV"

while read -r og; do
  # skip blanks/comments
  [[ -z "${og// }" ]] && continue
  [[ "$og" =~ ^# ]] && continue

  in="${OGDIR}/${og}.fasta"
  out="${OUTDIR}/${og}.fasta"

  if [ ! -s "$in" ]; then
    echo "[WARN] Missing: $in" | tee -a "$MISSING_LOG" >&2
    continue
  fi

  # Subset sequences by EXACT match to the species token after '|'
  # (headers look like: >g009075.m1|aspergillus_nidulans.final)
  awk '
    NR==FNR { sp[$1]=1; next }
    /^>/ {
      keep=0
      h=$0
      sub(/^>/,"",h)
      n=split(h,a,"|")
      if (n>=2) {
        s=a[2]
        sub(/[[:space:]].*$/,"",s)
        if (sp[s]) keep=1
      }
    }
    { if (keep) print }
  ' "$SP" "$in" > "$out"

  n=$(grep -c '^>' "$out" 2>/dev/null || true)

  echo -e "${og}\t${n}" >> "$COUNTS_TSV"

  if [ "${n}" -eq 0 ]; then
    echo "[INFO] ${og} has 0 sequences after filtering (check species naming vs headers): $out" | tee -a "$ZERO_LOG" >&2
  fi

done < "$LIST"

echo "[DONE] Subsetting complete."
echo "  Output FASTAs : $OUTDIR"
echo "  Counts TSV    : $COUNTS_TSV"
echo "  Missing OGs   : $MISSING_LOG"
echo "  Zero after    : $ZERO_LOG"
