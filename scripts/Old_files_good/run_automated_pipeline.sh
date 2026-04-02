#!/bin/bash
# run_automated_pipeline_all.sh - Complete automated pipeline for ALL 9 alignments

set -e  # Exit on any error

# Configuration
BASE_DIR="/mnt/bigdata/linuxhome/narhmadey/2_1.Y1000_267_phylogenetics/2_Plant_Path_563"
ALIGN_DIR="$BASE_DIR/data"
RAXML_SCRIPT="$BASE_DIR/scripts/run_raxml_advanced.sh"
SUBMIT_FILE="$BASE_DIR/scripts/raxml_jobs.sub"

# Create a directory for IQ-TREE logs
mkdir -p "$BASE_DIR/iqtree_logs"

# Loop through each alignment file
for ALIGNMENT in "$ALIGN_DIR"/OG*.trim.aln.fasta; do
    # Extract base filename without path and extension
    FILENAME=$(basename "$ALIGNMENT")
    BASENAME="${FILENAME%.fasta}"
    
    echo "========================================="
    echo "Processing: $FILENAME"
    echo "========================================="
    
    # Define prefix for this alignment's output
    IQTREE_PREFIX="$BASE_DIR/iqtree_logs/${BASENAME}_model"
    
    echo "=== STEP 1: Running IQ-TREE ModelFinder for $FILENAME ==="
    
    # Run IQ-TREE model selection
    iqtree -s "$ALIGNMENT" \
           -st AA \
           -m MF \
           -T AUTO \
           -pre "$IQTREE_PREFIX"
    
    echo "=== STEP 2: Extracting best model for $FILENAME ==="
    
    # Extract best model from IQ-TREE log
    BEST_MODEL=$(grep "Best-fit model" "${IQTREE_PREFIX}.log" | head -1 | sed 's/.*: \([^ ]*\) .*/\1/')
    
    if [ -z "$BEST_MODEL" ]; then
        echo "ERROR: Could not extract best model for $FILENAME"
        exit 1
    fi
    
    echo "Best model for $FILENAME: $BEST_MODEL"
    
    # Update RAxML script with the best model for this alignment
    # Note: This updates the script in place – we'll restore original later
    sed -i "s/^CANDIDATE_MODELS=\".*\"/CANDIDATE_MODELS=\"$BEST_MODEL\"/" "$RAXML_SCRIPT"
    
    echo "=== STEP 3: Submitting RAxML job for $FILENAME with model $BEST_MODEL ==="
    
    # Create a temporary input.txt with just this file
    echo "$FILENAME" > "$BASE_DIR/scripts/input.txt"
    
    cd "$BASE_DIR/scripts"
    condor_submit "$SUBMIT_FILE"
    
    # Wait a moment between submissions to avoid overwhelming the scheduler
    sleep 2
    
    echo "Job submitted for $FILENAME"
    echo ""
done

# Restore original RAxML script with placeholder model
sed -i 's/^CANDIDATE_MODELS=\".*\"/CANDIDATE_MODELS=\"LG+G4\"/' "$RAXML_SCRIPT"

echo "========================================="
echo "=== All 9 jobs submitted successfully ==="
echo "========================================="
echo "Check status with: condor_q"
echo "IQ-TREE logs are in: $BASE_DIR/iqtree_logs/"
echo "RAxML results will appear in: $BASE_DIR/results/"
echo "RAxML logs will appear in: $BASE_DIR/logs/"