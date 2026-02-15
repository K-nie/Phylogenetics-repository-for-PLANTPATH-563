#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# 01_qc_seqkit.sh
# Safe QC runner:
# - Creates missing folders with mkdir -p (never deletes anything)
# - Logs all commands
# - Runs seqkit stats on subset genomes + proteins (and CDS if present)
# How to run:chmod +x scripts/01_qc_seqkit.sh
#  ./scripts/01_qc_seqkit.sh
# Expected inputs (preferred):
#   data/processed/02_subset/genome/*
#   data/processed/02_subset/pep/*
# Optional:
#   data/processed/02_subset/cds/*
# ------------------------------------------------------------

# Project root = parent of scripts/
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DATA_DIR="${PROJECT_ROOT}/data"
PROC_DIR="${DATA_DIR}/processed"
SUBSET_DIR="${PROC_DIR}/02_subset"
QC_DIR="${PROC_DIR}/09_qc"
REPORT_DIR="${PROJECT_ROOT}/reports"

# Output files
GENOME_STATS="${QC_DIR}/assemblies/genome_seqkit_stats.tsv"
PEP_STATS="${QC_DIR}/proteins/pep_seqkit_stats.tsv"
CDS_STATS="${QC_DIR}/proteins/cds_seqkit_stats.tsv"

# Create folders (safe)
mkdir -p "${PROJECT_ROOT}/"{docs,scripts,results}
mkdir -p "${DATA_DIR}"/{raw,processed}
mkdir -p "${SUBSET_DIR}"/{genome,pep,cds}   # does NOT delete anything if already exists
mkdir -p "${QC_DIR}"/{assemblies,proteins,msa,logs,results}
mkdir -p "${REPORT_DIR}/logs"

# Log everything
LOG_FILE="${REPORT_DIR}/logs/qc_commands.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

date
echo "[INFO] Project root: ${PROJECT_ROOT}"
echo "[INFO] Subset dir:   ${SUBSET_DIR}"
echo "[INFO] QC dir:       ${QC_DIR}"
echo

# Helper: fail early if no files match a pattern
require_any() {
  local pattern="$1"
  shopt -s nullglob
  local files=( $pattern )
  shopt -u nullglob
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "[ERROR] No files matched: ${pattern}"
    echo "       Put your subset FASTA files into the expected folder."
    exit 1
  fi
}

# -----------------------------
# QC: genome assemblies
# -----------------------------
require_any "${SUBSET_DIR}/genome/*"
echo "[INFO] Running seqkit stats on genomes..."
seqkit stats -a "${SUBSET_DIR}/genome/"* > "${GENOME_STATS}"
echo "[INFO] Wrote: ${GENOME_STATS}"
head -n 5 "${GENOME_STATS}"
echo

# -----------------------------
# QC: proteins (pep)
# -----------------------------
require_any "${SUBSET_DIR}/pep/*"
echo "[INFO] Running seqkit stats on proteins (pep)..."
seqkit stats -a "${SUBSET_DIR}/pep/"* > "${PEP_STATS}"
echo "[INFO] Wrote: ${PEP_STATS}"
head -n 5 "${PEP_STATS}"
echo

# -----------------------------
# Optional: CDS
# -----------------------------
if compgen -G "${SUBSET_DIR}/cds/*" > /dev/null; then
  echo "[INFO] Running seqkit stats on CDS..."
  seqkit stats -a "${SUBSET_DIR}/cds/"* > "${CDS_STATS}"
  echo "[INFO] Wrote: ${CDS_STATS}"
  head -n 5 "${CDS_STATS}"
  echo
else
  echo "[INFO] No CDS files found in ${SUBSET_DIR}/cds/. Skipping CDS QC."
  echo
fi

echo "[DONE] QC complete."
date
