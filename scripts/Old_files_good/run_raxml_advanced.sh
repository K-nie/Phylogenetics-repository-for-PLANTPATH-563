#!/bin/bash
set -e
set -o pipefail

FILENAME="$1"

BASE_DIR="/mnt/bigdata/linuxhome/narhmadey/2_1.Y1000_267_phylogenetics/2_Plant_Path_563"
ALIGN_DIR="$BASE_DIR/data"                    # FIXED: alignments are in data/
RESULTS_DIR="$BASE_DIR/results"
LOGS_DIR="$BASE_DIR/logs"
INPUT_ALN="$ALIGN_DIR/$FILENAME"
OUTPUT_PREFIX="${FILENAME%.fasta}"            # remove .fasta extension

ENV_TAR="$BASE_DIR/phylo_env.tar.gz"

# List of candidate models (comma-separated). If only one, no model selection.
CANDIDATE_MODELS="LG"

# Number of bootstrap replicates
BOOTSTRAP_REPS=100

cd $_CONDOR_SCRATCH_DIR

echo "=== SCRATCH DIRECTORY CONTENTS (before unpack) ==="
ls -la
echo "=== END ==="

if [ ! -f "$ENV_TAR" ]; then
    echo "ERROR: Environment tarball $ENV_TAR not found"
    exit 1
fi

cp "$ENV_TAR" ./
tar -xzf phylo_env.tar.gz
rm -f phylo_env.tar.gz

export PATH="$(pwd)/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd)/lib:$(pwd)/lib64:$LD_LIBRARY_PATH"

BINARY="$(pwd)/bin/raxml-ng"
if [ ! -x "$BINARY" ]; then
    echo "ERROR: Binary $BINARY not found or not executable"
    exit 1
fi

echo "=== RAXML BINARY INFO ==="
ldd "$BINARY" || true
echo "=== END ==="

# ------------------------------------------------------------
# Determine if we have multiple models or just one
# ------------------------------------------------------------
if [[ "$CANDIDATE_MODELS" == *","* ]]; then
    # Multiple models: run model selection
    echo "=== STEP 1: Model selection with candidate models: $CANDIDATE_MODELS ==="
    "$BINARY" --msa "$INPUT_ALN" \
              --data-type AA \
              --model "$CANDIDATE_MODELS" \
              --prefix "model_sel" \
              --threads "${OMP_NUM_THREADS:-1}" \
              --seed 12345 \
              --opt-model BIC \
              --redo \
              2>&1 | tee model_selection.log

    if [ ! -f "model_sel.raxml.bestModel" ]; then
        echo "ERROR: Model selection did not produce model_sel.raxml.bestModel"
        exit 1
    fi

    BEST_MODEL=$(cat model_sel.raxml.bestModel | tr -d '\n')
    echo "=== Best model selected: $BEST_MODEL ==="
else
    # Single model: no selection needed
    BEST_MODEL="$CANDIDATE_MODELS"
    echo "=== Using fixed model: $BEST_MODEL ==="
fi

# ------------------------------------------------------------
# STEP 2: Full tree search + bootstrap with best model
# ------------------------------------------------------------
echo "=== STEP 2: Full tree search and bootstrap with model $BEST_MODEL ==="
"$BINARY" --msa "$INPUT_ALN" \
          --data-type AA \
          --model "$BEST_MODEL" \
          --prefix "$OUTPUT_PREFIX" \
          --threads "${OMP_NUM_THREADS:-1}" \
          --seed 12345 \
          --tree pars{10},rand{10} \
          --bs-trees $BOOTSTRAP_REPS \
          --all \
          2>&1 | tee "$OUTPUT_PREFIX.raxml.log"

# ------------------------------------------------------------
# Copy results to shared directories
# ------------------------------------------------------------
echo "=== Copying results to $RESULTS_DIR and $LOGS_DIR ==="
mkdir -p "$RESULTS_DIR" "$LOGS_DIR"
cp "$OUTPUT_PREFIX.raxml.bestTree" "$RESULTS_DIR/" 2>/dev/null || echo "Warning: bestTree not found"
cp "$OUTPUT_PREFIX.raxml.bestModel" "$RESULTS_DIR/" 2>/dev/null || echo "Warning: bestModel not found"
cp "$OUTPUT_PREFIX.raxml.bootstraps" "$RESULTS_DIR/" 2>/dev/null || echo "Warning: bootstraps not found"
cp "$OUTPUT_PREFIX.raxml.support" "$RESULTS_DIR/" 2>/dev/null || echo "Warning: support not found"
cp "$OUTPUT_PREFIX.raxml.log" "$LOGS_DIR/" 2>/dev/null || echo "Warning: log not found"

echo "=== All done ==="