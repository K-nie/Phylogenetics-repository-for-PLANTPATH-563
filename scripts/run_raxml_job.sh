#!/bin/bash
# run_raxml_job.sh
# Condor job script: runs RAxML-NG tree search for a single orthogroup.
# Called by raxml_jobs.sub with two arguments:
#   $1 = alignment filename (e.g. OG0000007.trim.aln.fasta)
#   $2 = best-fit model string (e.g. Q.YEAST+F+R10)
#
# The model is passed as an argument — NOT read from a shared file —
# so all jobs can run simultaneously without race conditions.
#Author:Benjamin Narh-Madey

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
FILENAME="$1"
BEST_MODEL="$2"

if [ -z "$FILENAME" ] || [ -z "$BEST_MODEL" ]; then
    echo "ERROR: Usage: run_raxml_job.sh <alignment_filename> <model>"
    exit 1
fi

# ── Paths (absolute — Condor execute nodes may have different CWD) ─────────────
BASE_DIR="/home/glbrc.org/narhmadey/2.Y1000_267_phylogenetics/2_Plant_Path_563"
ALIGN_DIR="$BASE_DIR/data"
RESULTS_DIR="$BASE_DIR/results"
LOGS_DIR="$BASE_DIR/logs"
ENV_TAR="$BASE_DIR/phylo_env.tar.gz"

INPUT_ALN="$ALIGN_DIR/$FILENAME"
OUTPUT_PREFIX="${FILENAME%.fasta}"      # e.g. OG0000007.trim.aln

echo "========================================="
echo "Job started:   $(date)"
echo "Host:          $(hostname)"
echo "Alignment:     $FILENAME"
echo "Model:         $BEST_MODEL"
echo "OMP threads:   ${OMP_NUM_THREADS:-1}"
echo "Scratch dir:   $_CONDOR_SCRATCH_DIR"
echo "========================================="

# ── Validate input ────────────────────────────────────────────────────────────
if [ ! -f "$INPUT_ALN" ]; then
    echo "ERROR: Alignment not found: $INPUT_ALN"
    exit 1
fi

# ── Set up environment in scratch ─────────────────────────────────────────────
cd "$_CONDOR_SCRATCH_DIR"

echo "=== Unpacking environment ==="
if [ ! -f "$ENV_TAR" ]; then
    echo "ERROR: Environment tarball not found: $ENV_TAR"
    exit 1
fi
cp "$ENV_TAR" ./
tar -xzf phylo_env.tar.gz
rm -f phylo_env.tar.gz

export PATH="$(pwd)/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd)/lib:$(pwd)/lib64:${LD_LIBRARY_PATH:-}"

# Prefer the pre-built raxml-ng from bin/ in the project if env unpack lacks it
RAXML_BIN="$(pwd)/bin/raxml-ng"
if [ ! -x "$RAXML_BIN" ]; then
    # Fall back to project bin/
    RAXML_BIN="$BASE_DIR/bin/raxml-ng"
fi
if [ ! -x "$RAXML_BIN" ]; then
    echo "ERROR: raxml-ng not found in scratch or $BASE_DIR/bin/"
    exit 1
fi

echo "=== RAxML-NG binary: $RAXML_BIN ==="
"$RAXML_BIN" --version

# ── Step 1: Evaluate parsimony starting trees ─────────────────────────────────
echo ""
echo "=== STEP 1: Tree search (model: $BEST_MODEL) ==="
"$RAXML_BIN" \
    --all \
    --msa "$INPUT_ALN" \
    --data-type AA \
    --model "$BEST_MODEL" \
    --prefix "$OUTPUT_PREFIX" \
    --threads "${OMP_NUM_THREADS:-1}" \
    --seed 12345 \
    --tree pars{10},rand{10} \
    --bs-trees 1000 \
    --redo \
    2>&1 | tee "${OUTPUT_PREFIX}.raxml.screen.log"

# ── Copy results to shared filesystem ────────────────────────────────────────
echo ""
echo "=== Copying results to $RESULTS_DIR ==="
mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

# Core result files
for EXT in bestTree bestModel support bootstraps; do
    SRC="${OUTPUT_PREFIX}.raxml.${EXT}"
    if [ -f "$SRC" ]; then
        cp "$SRC" "$RESULTS_DIR/"
        echo "  Copied: $SRC"
    else
        echo "  Warning: $SRC not found (may be normal for some run modes)"
    fi
done

# Logs
cp "${OUTPUT_PREFIX}.raxml.screen.log" "$LOGS_DIR/" 2>/dev/null || true
cp "${OUTPUT_PREFIX}.raxml.log"        "$LOGS_DIR/" 2>/dev/null || true

echo ""
echo "========================================="
echo "Job finished:  $(date)"
echo "Results in:    $RESULTS_DIR"
echo "========================================="