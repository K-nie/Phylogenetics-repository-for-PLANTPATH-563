#!/usr/bin/env bash
set -euo pipefail

# Paths
LIST="data/processed/species_267.txt"
TARBALL_DIR="data/processed/01_extract/y1000p_cds_files/y1000p_cds_files"
OUT_DIR="data/processed/02_subset"
OUT_CDS="${OUT_DIR}/cds"
LOG_DIR="data/processed/09_qc/logs"
LOG_FILE="${LOG_DIR}/subset_267_extract.log"
FOUND_LIST="${LOG_DIR}/subset_267_found_tarballs.txt"
MISSING_LIST="${LOG_DIR}/subset_267_missing.txt"

mkdir -p "${OUT_CDS}" "${LOG_DIR}"

# OPTIONAL but recommended: clear previous subset outputs (not raw data)
rm -f "${OUT_CDS}/"* 2>/dev/null || true

# Log everything
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[INFO] Using species list: ${LIST}"
echo "[INFO] Searching tarballs in: ${TARBALL_DIR}"
echo "[INFO] Output CDS folder: ${OUT_CDS}"
echo

# Basic checks
if [[ ! -f "${LIST}" ]]; then
  echo "[ERROR] Cannot find species list file: ${LIST}"
  exit 1
fi
if [[ ! -d "${TARBALL_DIR}" ]]; then
  echo "[ERROR] Cannot find tarball directory: ${TARBALL_DIR}"
  exit 1
fi

# Clean previous "found/missing" reports (not your data)
: > "${FOUND_LIST}"
: > "${MISSING_LIST}"

echo "[INFO] Total tarballs available: $(find "${TARBALL_DIR}" -maxdepth 1 -type f -name "*.final.cds.tar.gz" | wc -l | tr -d ' ')"
echo

# For each line in species list, find matching tarball(s) by substring
while IFS= read -r raw || [[ -n "${raw}" ]]; do
  s="$(echo "${raw}" | sed 's/#.*//; s/^[[:space:]]*//; s/[[:space:]]*$//')"
  [[ -z "${s}" ]] && continue

  # Find matches (substring search against filenames)
  matches="$(find "${TARBALL_DIR}" -maxdepth 1 -type f -name "*.final.cds.tar.gz" -print | \
            awk -v pat="${s}" '
              BEGIN{ found=0 }
              {
                n=$0
                sub(/^.*\//,"",n)   # basename
                if (index(n, pat) > 0) { print $0; found=1 }
              }
              END{ if (found==0) exit 1 }
            ' 2>/dev/null || true)"

  if [[ -z "${matches}" ]]; then
    echo "${s}" >> "${MISSING_LIST}"
  else
    printf "%s\n" "${matches}" >> "${FOUND_LIST}"
  fi
done < "${LIST}"

# Deduplicate found list
sort -u "${FOUND_LIST}" -o "${FOUND_LIST}"

FOUND_N="$(wc -l < "${FOUND_LIST}" | tr -d ' ')"
MISS_N="$(wc -l < "${MISSING_LIST}" | tr -d ' ')"

echo "[INFO] Matched tarballs: ${FOUND_N}"
echo "[INFO] Missing entries:  ${MISS_N}"
echo "[INFO] Found list:       ${FOUND_LIST}"
echo "[INFO] Missing list:     ${MISSING_LIST}"
echo

if [[ "${FOUND_N}" -eq 0 ]]; then
  echo "[ERROR] No matches found. Your species_267.txt strings may not match filenames."
  echo "        Inspect one filename in ${TARBALL_DIR} and one entry in ${LIST} and align them."
  exit 1
fi

# Extract only the matched tarballs into OUT_CDS
echo "[INFO] Extracting matched tarballs into: ${OUT_CDS}"
i=0
while IFS= read -r tarf; do
  i=$((i+1))
  echo "[INFO] (${i}/${FOUND_N}) Extract: $(basename "$tarf")"
  tar -xzf "$tarf" -C "${OUT_CDS}"
done < "${FOUND_LIST}"

echo
echo "[DONE] Subset extraction complete."
echo "[INFO] Example output files:"
ls -lah "${OUT_CDS}" | head -n 20
