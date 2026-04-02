#!/bin/bash
# collect_models.sh
# Runs on the submit node AFTER all ModelFinder Condor jobs complete.
# Scans iqtree_logs/ for best-fit models and writes jobs_with_models.txt
# which is then used by raxml_jobs.sub and beast_jobs.sub.
#
# USAGE:
#   bash scripts/collect_models.sh
#
# PREREQUISITES:
#   All modelfinder_jobs.sub jobs must be complete (condor_q shows empty)
#
# OUTPUT:
#   scripts/jobs_with_models.txt  — tab-separated: filename<TAB>best_model
#Author:Benjamin Narh-Madey

set -euo pipefail

BASE_DIR="/home/glbrc.org/narhmadey/2.Y1000_267_phylogenetics/2_Plant_Path_563"
ALIGN_DIR="$BASE_DIR/data"
IQTREE_LOGS="$BASE_DIR/iqtree_logs"
JOBS_FILE="$BASE_DIR/scripts/jobs_with_models.txt"

echo "========================================="
echo "Collecting ModelFinder results"
echo "Started: $(date)"
echo "========================================="

# Clear jobs file
> "$JOBS_FILE"

TOTAL=0
FAILED=0
FAILED_LIST=()

for ALIGNMENT in "$ALIGN_DIR"/OG*.trim.aln.fasta; do
    FILENAME=$(basename "$ALIGNMENT")
    BASENAME="${FILENAME%.fasta}"
    SCREEN_LOG="$IQTREE_LOGS/${BASENAME}_model.screen.log"
    IQTREE_LOG="$IQTREE_LOGS/${BASENAME}_model.log"

    echo ""
    echo "--- $FILENAME ---"

    # Try screen.log first, fall back to .log
    LOG_TO_USE=""
    if [ -f "$SCREEN_LOG" ] && grep -q "Best-fit model:" "$SCREEN_LOG"; then
        LOG_TO_USE="$SCREEN_LOG"
    elif [ -f "$IQTREE_LOG" ] && grep -q "Best-fit model:" "$IQTREE_LOG"; then
        LOG_TO_USE="$IQTREE_LOG"
    fi

    if [ -z "$LOG_TO_USE" ]; then
        echo "  ERROR: No completed ModelFinder log found"
        echo "         Expected: $SCREEN_LOG"
        echo "         Check if the Condor job completed successfully"
        FAILED=$((FAILED + 1))
        FAILED_LIST+=("$FILENAME")
        continue
    fi

    # Extract best model
    BEST_MODEL=$(grep "Best-fit model:" "$LOG_TO_USE" \
                 | head -1 \
                 | sed 's/.*Best-fit model: \([^ ]*\) .*/\1/')

    if [ -z "$BEST_MODEL" ]; then
        echo "  ERROR: Could not extract best model from $LOG_TO_USE"
        FAILED=$((FAILED + 1))
        FAILED_LIST+=("$FILENAME")
        continue
    fi

    echo "  Best model: $BEST_MODEL"
    printf "%s\t%s\n" "$FILENAME" "$BEST_MODEL" >> "$JOBS_FILE"
    TOTAL=$((TOTAL + 1))
done

echo ""
echo "========================================="
echo "Collection complete: $(date)"
echo "  Successful: $TOTAL"
echo "  Failed:     $FAILED"
echo "  Jobs file:  $JOBS_FILE"
echo "========================================="

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "WARNING: The following orthogroups failed:"
    for f in "${FAILED_LIST[@]}"; do
        echo "  - $f"
    done
    echo ""
    echo "To resubmit failed jobs:"
    echo "  1. Check Condor logs in logs/mf_*.error for error details"
    echo "  2. Recreate modelfinder_jobs.txt with only failed OGs:"
    for f in "${FAILED_LIST[@]}"; do
        echo "     echo $f >> scripts/modelfinder_jobs_retry.txt"
    done
    echo "  3. condor_submit scripts/modelfinder_jobs.sub"
    exit 1
fi

echo ""
echo "All models collected. Review jobs_with_models.txt:"
echo "  cat $JOBS_FILE"
echo ""
echo "Then submit Phase 2 jobs:"
echo "  condor_submit $BASE_DIR/scripts/raxml_jobs.sub"
echo "  python3 $BASE_DIR/scripts/generate_beast_xml.py"
echo "  ls beast_xml/*.xml | xargs -n1 basename | sed 's/\\.xml//' > $BASE_DIR/scripts/beast_jobs.txt"
echo "  condor_submit $BASE_DIR/scripts/beast_jobs.sub"