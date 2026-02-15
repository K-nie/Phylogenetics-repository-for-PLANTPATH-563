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

# Log everything
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[INFO] Using species list: ${LIST}"
echo "[INFO] Searching tarballs in: ${TARBALL_DIR}"
echo "[INFO] Output CDS folder: ${OUT_CDS}"
echo

# Basic checks
if [[ ! -f "${LIST}" ]]; then
  echo "[ERROR] Cannot find ${LIST} in the project root."
  exit 1
fi
if [[ ! -d "${TARBALL_DIR}" ]]; then
  echo "[ERROR] Cannot find tarball directory: ${TARBALL_DIR}"
  exit 1
fi

# Clean previous "found/missing" reports (not your data)
: > "${FOUND_LIST}"
: > "${MISSING_LIST}"

# Build a filename inventory once (fast)
# Example tarball names:
# yHMPu5000026055_diutina_siamensis_170912.final.cds.tar.gz
mapfile -t ALL_TARS < <(find "${TARBALL_DIR}" -maxdepth 1 -type f -name "*.final.cds.tar.gz" -print)

echo "[INFO] Total tarballs available: ${#ALL_TARS[@]}"
echo

# For each line in species_267.txt, try to find matching tarball(s)
# We match by substring, so your list can contain:
#  - full prefix (yHMPu...._species_strain)
#  - species name fragment (diutina_siamensis)
#  - exact tarball stem without extension
while IFS= read -r raw || [[ -n "${raw}" ]]; do
  s="$(echo "${raw}" | sed 's/#.*//; s/^[[:space:]]*//; s/[[:space:]]*$//')"
  [[ -z "${s}" ]] && continue

  # Find tarballs whose filename contains the species token
  matches=()
  for f in "${ALL_TARS[@]}"; do
    base="$(basename "$f")"
    if [[ "$base" == *"$s"* ]]; then
      matches+=("$f")
    fi
  done

  if [[ ${#matches[@]} -eq 0 ]]; then
    echo "$s" >> "${MISSING_LIST}"
  else
    # If multiple match, we take all (you can tighten later if needed)
    for m in "${matches[@]}"; do
      echo "$m" >> "${FOUND_LIST}"
    done
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
  ((i+=1))
  echo "[INFO] (${i}/${FOUND_N}) Extract: $(basename "$tarf")"
  tar -xzf "$tarf" -C "${OUT_CDS}"
done < "${FOUND_LIST}"

echo
echo "[DONE] Subset extraction complete."
echo "[INFO] Example output files:"
ls -lah "${OUT_CDS}" | head -n 20
