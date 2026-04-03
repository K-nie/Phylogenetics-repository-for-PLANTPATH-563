#!/bin/bash
# run_beast_job.sh
# Condor job script: runs BEAST2 Bayesian inference for a single orthogroup.
# Called by beast_jobs.sub with one argument:
#   $1 = orthogroup name (e.g. OG0000007)
#
# Settings: Relaxed lognormal clock + Birth-Death tree prior + 50M MCMC
#Author:Benjamin Narh-Madey


set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────
OG_NAME="$1"

if [ -z "$OG_NAME" ]; then
    echo "ERROR: Usage: run_beast_job.sh <og_name>"
    echo "       e.g.   run_beast_job.sh OG0000007"
    exit 1
fi

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR="/home/glbrc.org/narhmadey/2.Y1000_267_phylogenetics/2_Plant_Path_563"
XML_DIR="$BASE_DIR/beast_xml"
RESULTS_DIR="$BASE_DIR/results/beast"
LOGS_DIR="$BASE_DIR/logs"
ENV_TAR="$BASE_DIR/phylo_env.tar.gz"
XML_FILE="$XML_DIR/${OG_NAME}.xml"

echo "========================================="
echo "BEAST2 Job started:  $(date)"
echo "Host:                $(hostname)"
echo "Orthogroup:          $OG_NAME"
echo "XML config:          $XML_FILE"
echo "Scratch dir:         $_CONDOR_SCRATCH_DIR"
echo "========================================="

# ── Validate inputs ───────────────────────────────────────────────────────────
if [ ! -f "$XML_FILE" ]; then
    echo "ERROR: XML file not found: $XML_FILE"
    echo "       Run generate_beast_xml.py first"
    exit 1
fi

# ── Set up environment in scratch ─────────────────────────────────────────────
cd "$_CONDOR_SCRATCH_DIR"

echo "=== Unpacking environment ==="
cp "$ENV_TAR" ./
tar -xzf phylo_env.tar.gz
rm -f phylo_env.tar.gz

export PATH="$(pwd)/bin:$PATH"
export LD_LIBRARY_PATH="$(pwd)/lib:$(pwd)/lib64:${LD_LIBRARY_PATH:-}"

# BEAST2 needs Java — use the one from the unpacked env
export JAVA_HOME="$(pwd)"
BEAST_BIN="$(pwd)/bin/beast"

if [ ! -x "$BEAST_BIN" ]; then
    echo "ERROR: beast binary not found in unpacked environment"
    echo "       Expected: $BEAST_BIN"
    exit 1
fi

echo "=== BEAST2 binary: $BEAST_BIN ==="
"$BEAST_BIN" -version

# ── Copy XML to scratch ───────────────────────────────────────────────────────
cp "$XML_FILE" ./

# ── Run BEAST2 ────────────────────────────────────────────────────────────────
echo ""
echo "=== Running BEAST2 (50M MCMC chain) ==="
echo "=== Relaxed lognormal clock + Birth-Death prior ==="
echo "=== Seed: ${RANDOM_SEED:=$RANDOM} ==="
echo ""

# We use -seed $RANDOM so Run 1 and Run 2 start at different points
"$BEAST_BIN" \
    -beagle_off \
    -threads "${OMP_NUM_THREADS:-4}" \
    -overwrite \
    -seed "$RANDOM_SEED" \
    "${OG_NAME}.xml" \
    2>&1 | tee "${OG_NAME}.beast.screen.log"

# ── Copy results to shared filesystem ────────────────────────────────────────
echo ""
echo "=== Copying results to $RESULTS_DIR ==="
mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

# Core BEAST2 output files
# generate_beast_xml.py names files as {OG_NAME}_{run_id}.trace.log / .trees
# The run_id is embedded in the XML passed to this job — match by glob
for EXT in trace.log trees; do
    for SRC in "${OG_NAME}"_run*.${EXT}; do
        if [ -f "$SRC" ]; then
            cp "$SRC" "$RESULTS_DIR/"
            echo "  Copied: $SRC"
        fi
    done
    # Warn only if nothing matched at all
    found=$(ls "${OG_NAME}"_run*.${EXT} 2>/dev/null | wc -l)
    if [ "$found" -eq 0 ]; then
        echo "  Warning: no ${OG_NAME}_run*.${EXT} files found"
    fi
done

# Screen log
cp "${OG_NAME}.beast.screen.log" "$LOGS_DIR/" 2>/dev/null || true

echo ""
echo "========================================="
echo "BEAST2 Job finished: $(date)"
echo "Results in:          $RESULTS_DIR"
echo "========================================="
echo ""
echo "Next steps after all jobs complete:"
echo "  1. Check convergence in Tracer (ESS > 200 for all parameters)"
echo "  2. Summarize posterior trees with TreeAnnotator:"
echo "     treeannotator -burnin 10 ${OG_NAME}.trees ${OG_NAME}.beast.tree"