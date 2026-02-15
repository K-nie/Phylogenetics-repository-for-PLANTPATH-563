#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Trim MAFFT alignments with trimAl (robust + logged)
# -----------------------------

# Avoid "unbound variable" crashes from shells/plugins referencing HISTTIMEFORMAT
: "${HISTTIMEFORMAT:=}"

IN="data/processed/03_ml_orthogroups/06_alignments"
OUT="data/processed/03_ml_orthogroups/07_trimmed_alignments"
LOG="data/processed/03_ml_orthogroups/logs/trimal_run.log"

mkdir -p "$OUT" "$(dirname "$LOG")"
: > "$LOG"

# Dependency check
command -v trimal >/dev/null 2>&1 || { echo "[ERROR] trimal not found in PATH" | tee -a "$LOG" >&2; exit 1; }

# Fail-fast if no inputs match
shopt -s nullglob
files=("$IN"/*.aln.fasta)
shopt -u nullglob

echo "[INFO] Inputs matched: ${#files[@]} files" | tee -a "$LOG"
if [ ${#files[@]} -eq 0 ]; then
  echo "[ERROR] No MAFFT alignments found: $IN/*.aln.fasta" | tee -a "$LOG" >&2
  exit 1
fi

for aln in "${files[@]}"; do
  og=$(basename "$aln" .aln.fasta)
  out="${OUT}/${og}.trim.aln.fasta"
  tmp="${out}.tmp"

  echo "[TRIMAL] $og" | tee -a "$LOG"

  # Run trimAl (log everything), write to temp then atomically move into place
  if command -v /usr/bin/time >/dev/null 2>&1; then
    /usr/bin/time -l trimal -in "$aln" -out "$tmp" -automated1 >>"$LOG" 2>&1
  else
    trimal -in "$aln" -out "$tmp" -automated1 >>"$LOG" 2>&1
  fi
  mv -f "$tmp" "$out"

  # Sanity check
  if [ ! -s "$out" ] || ! grep -q '^>' "$out"; then
    echo "[ERROR] Output empty for $og: $out" | tee -a "$LOG" >&2
    exit 1
  fi
done

echo "[DONE] TrimAl complete: wrote ${#files[@]} files to $OUT" | tee -a "$LOG"
