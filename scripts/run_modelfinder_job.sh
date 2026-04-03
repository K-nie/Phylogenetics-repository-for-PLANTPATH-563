#!/bin/bash
# Author:Benjamin Narh-Madey
# run_modelfinder_job.sh
# Condor job script: runs IQ-TREE ModelFinder for a single orthogroup.
#
set -euo pipefail

FILENAME="$1"
if [ -z "$FILENAME" ]; then
    echo "ERROR: Usage: run_modelfinder_job.sh <alignment_filename>"
    exit 1
fi

BASE_DIR="/home/glbrc.org/narhmadey/2.Y1000_267_phylogenetics/2_Plant_Path_563"
ALIGN_DIR="$BASE_DIR/data"
IQTREE_LOGS="$BASE_DIR/iqtree_logs"
ENV_TAR="$BASE_DIR/phylo_env.tar.gz"
INPUT_ALN="$ALIGN_DIR/$FILENAME"
BASENAME="${FILENAME%.fasta}"
OUTPUT_PREFIX="${BASENAME}_model"

echo "========================================="
echo "ModelFinder Job started: $(date)"
echo "Host:                    $(hostname)"
echo "Alignment:               $FILENAME"
echo "Scratch dir:             $_CONDOR_SCRATCH_DIR"
echo "========================================="

if [ ! -f "$INPUT_ALN" ]; then
    echo "ERROR: Alignment not found: $INPUT_ALN"
    exit 1
fi

cd "$_CONDOR_SCRATCH_DIR"

echo "=== Unpacking environment ==="
cp "$ENV_TAR" ./
tar -xzf phylo_env.tar.gz
rm -f phylo_env.tar.gz

export PATH="$(pwd)/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd)/lib:$(pwd)/lib64:${LD_LIBRARY_PATH:-}"

IQTREE_BIN="$(pwd)/bin/iqtree"
if [ ! -x "$IQTREE_BIN" ]; then
    IQTREE_BIN="$(which iqtree 2>/dev/null || echo '')"
fi
if [ -z "$IQTREE_BIN" ] || [ ! -x "$IQTREE_BIN" ]; then
    echo "ERROR: iqtree binary not found"
    exit 1
fi

echo "=== Running IQ-TREE ModelFinder ==="
"$IQTREE_BIN" \
    -s "$INPUT_ALN" \
    -st AA \
    -m MFP \
    -B 1000 \
    -T AUTO \
    --threads-max "${OMP_NUM_THREADS:-4}" \
    -pre "$OUTPUT_PREFIX" \
    --redo \
    2>&1 | tee "${OUTPUT_PREFIX}.screen.log"
    
if ! grep -q "Best-fit model:" "${OUTPUT_PREFIX}.screen.log"; then
    echo "ERROR: ModelFinder did not produce a best-fit model"
    exit 1
fi

BEST_MODEL=$(grep "Best-fit model:" "${OUTPUT_PREFIX}.screen.log" \
             | head -1 \
             | sed 's/.*Best-fit model: \([^ ]*\) .*/\1/')

echo ""
echo "========================================="
echo "Best-fit model: $BEST_MODEL"
echo "========================================="

mkdir -p "$IQTREE_LOGS"
for EXT in log screen.log iqtree model.gz treefile uniqueseq.phy; do
    SRC="${OUTPUT_PREFIX}.${EXT}"
    if [ -f "$SRC" ]; then
        cp "$SRC" "$IQTREE_LOGS/"
        echo "  Copied: $SRC"
    fi
done

echo ""
echo "========================================="
echo "ModelFinder Job finished: $(date)"
echo "========================================="
