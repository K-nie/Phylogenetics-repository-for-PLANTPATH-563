#!/usr/bin/env bash
set -eo pipefail   # no 'u' → avoids unbound variable errors

# ──────────────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────────────

LIST="data/processed/species_267.txt"
RAW_TARBALL="data/raw/y1000p_gff3_files.tar.gz"
OUT_DIR="data/processed/10_genomeviz/05_subset/gff3"
LOG_DIR="data/processed/10_1_qc/logs"

LOG_FILE="${LOG_DIR}/subset_267_gff3.log"
FOUND_LIST="${LOG_DIR}/subset_267_gff3_found.txt"
MISSING_LIST="${LOG_DIR}/subset_267_gff3_missing.txt"
UNMATCHED_DIR="${OUT_DIR}/../unmatched_gff3_$(date +%Y%m%d_%H%M%S)"

mkdir -p "${OUT_DIR}" "${LOG_DIR}" "${UNMATCHED_DIR}"

# ──────────────────────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────────────────────

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "[START] $(date '+%Y-%m-%d %H:%M:%S')"
echo "[CONFIG] Species list : ${LIST}"
echo "[CONFIG] Tarball      : ${RAW_TARBALL}"
echo "[CONFIG] Output dir   : ${OUT_DIR}"
echo

# ──────────────────────────────────────────────────────────────────────────────
# Basic checks
# ──────────────────────────────────────────────────────────────────────────────

[[ -f "${LIST}" ]]         || { echo "[ERROR] Cannot find ${LIST}"; exit 1; }
[[ -f "${RAW_TARBALL}" ]]  || { echo "[ERROR] Cannot find ${RAW_TARBALL}"; exit 1; }

# ──────────────────────────────────────────────────────────────────────────────
# Extract once if needed
# ──────────────────────────────────────────────────────────────────────────────

if [[ ! -d "${OUT_DIR}/y1000p_gff3_files" || -z "$(ls -A "${OUT_DIR}/y1000p_gff3_files" 2>/dev/null)" ]]; then
    echo "[EXTRACTION] Extracting main tarball..."
    tar -xzf "${RAW_TARBALL}" -C "${OUT_DIR}"
    echo "[EXTRACTION] Done."
else
    echo "[INFO] Files already extracted → skipping"
fi
echo

# ──────────────────────────────────────────────────────────────────────────────
# Inventory of GFF3 files
# ──────────────────────────────────────────────────────────────────────────────

ALL_GFF3=()
while IFS= read -r f; do
    ALL_GFF3+=("$f")
done < <(find "${OUT_DIR}" -type f -name "*.final.gff3" -print | sort)

echo "[INFO] Total GFF3 files available: ${#ALL_GFF3[@]}"
(( ${#ALL_GFF3[@]} == 0 )) && { echo "[ERROR] No *.final.gff3 files found"; exit 1; }
echo

# ──────────────────────────────────────────────────────────────────────────────
# Match loop (no associative array)
# ──────────────────────────────────────────────────────────────────────────────

: > "${FOUND_LIST}"
: > "${MISSING_LIST}"

assigned_paths=()   # regular array instead of assoc array

while IFS= read -r raw || [[ -n "${raw}" ]]; do
    s="$(echo "${raw}" | sed 's/#.*//; s/^[[:space:]]*//; s/[[:space:]]*$//')"
    [[ -z "${s}" ]] && continue

    s_norm="$(echo "${s}" | tr '[:upper:]' '[:lower:]' | tr ' ' '_')"

    matches=()
    for f in "${ALL_GFF3[@]}"; do
        # Skip if already assigned
        if printf '%s\n' "${assigned_paths[@]}" | grep -Fxq "$f"; then
            continue
        fi

        base="$(basename "$f" .final.gff3)"
        base_norm="${base,,}"

        # Strip common prefixes
        cleaned="${base_norm}"
        cleaned="${cleaned#y[a-z]*[0-9]*_}"
        cleaned="${cleaned#yh[a-z0-9]*_}"

        # Match if species is prefix, exact, or bounded by _
        if [[ "${cleaned}" == "${s_norm}"* ]] ||
           [[ "${cleaned}" == "${s_norm}" ]] ||
           [[ "${cleaned}" =~ ^${s_norm}_ ]] ||
           [[ "${cleaned}" =~ _${s_norm}_ ]] ||
           [[ "${cleaned}" =~ _${s_norm}$ ]]; then

            matches+=("$f")
        fi
    done

    if [[ ${#matches[@]} -eq 0 ]]; then
        echo "$s" >> "${MISSING_LIST}"
    else
        first_match="${matches[0]}"
        echo "$first_match" >> "${FOUND_LIST}"
        assigned_paths+=("$first_match")

        if [[ ${#matches[@]} -gt 1 ]]; then
            echo "[WARN] $s matched ${#matches[@]} files → kept first one only" >> "${LOG_FILE}"
        fi
    fi
done < "${LIST}"

sort -u "${FOUND_LIST}" -o "${FOUND_LIST}"

FOUND_N=$(wc -l < "${FOUND_LIST}")
MISS_N=$(wc -l < "${MISSING_LIST}")

echo "[SUMMARY]"
echo "  Matched files    : ${FOUND_N}"
echo "  Missing species  : ${MISS_N}"
echo "  Found list       : ${FOUND_LIST}"
echo "  Missing list     : ${MISSING_LIST}"
echo

# ──────────────────────────────────────────────────────────────────────────────
# Cleanup unmatched
# ──────────────────────────────────────────────────────────────────────────────

if (( FOUND_N < ${#ALL_GFF3[@]} )); then
    echo "[CLEANUP] Moving unmatched files → ${UNMATCHED_DIR}"
    moved=0
    for f in "${ALL_GFF3[@]}"; do
        if ! grep -Fxq "$f" "${FOUND_LIST}"; then
            mv "$f" "${UNMATCHED_DIR}/"
            ((moved++))
        fi
    done
    echo "  → Moved ${moved} files"
else
    echo "[CLEANUP] All files matched → no cleanup"
fi

echo
echo "[DONE] Subset complete."
echo "[INFO] Kept files in: ${OUT_DIR}"
ls -lh "${OUT_DIR}" | head -n 15 || true
echo
echo "Log: ${LOG_FILE}"
echo "Finished: $(date '+%Y-%m-%d %H:%M:%S')"