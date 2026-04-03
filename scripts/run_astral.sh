#!/bin/bash
# Author:Benjamin Narh-Madey
# run_astral.sh
# Combines all RAxML-NG gene trees into a single multi-species coalescent tree.
# Author: Benjamin Narh-Madey

set -euo pipefail

# ── 0. Initialize Environment ─────────────────────────────────────────────────
# It is highly recommended to activate your conda environment inside the script 
# so the compute nodes know where the tools are.
source ~/.bashrc
conda activate phylo_env

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR="/home/glbrc.org/narhmadey/2.Y1000_267_phylogenetics/2_Plant_Path_563"
RAXML_RESULTS="$BASE_DIR/results"
ASTRAL_OUT_DIR="$BASE_DIR/results/astral"
LOGS_DIR="$BASE_DIR/logs"

echo "========================================="
echo "ASTRAL Species Tree Estimation Started"
echo "Time: $(date)"
echo "========================================="

# ── 1. Create Output Directories ──────────────────────────────────────────────
mkdir -p "$ASTRAL_OUT_DIR" "$LOGS_DIR"

# ── 2. Concatenate All Gene Trees ─────────────────────────────────────────────
# ASTRAL requires a single file containing all the newick trees, one per line.
COMBINED_TREES="$ASTRAL_OUT_DIR/all_raxml_gene_trees.tre"

echo "Gathering RAxML best trees..."
# We use .raxml.bestTree because ASTRAL expects unrooted trees
cat "$RAXML_RESULTS"/*.raxml.bestTree > "$COMBINED_TREES"

TREE_COUNT=$(wc -l < "$COMBINED_TREES")
echo "Successfully combined $TREE_COUNT gene trees into $COMBINED_TREES"

# ── 3. Run ASTRAL ─────────────────────────────────────────────────────────────
SPECIES_TREE="$ASTRAL_OUT_DIR/astral_species_tree.tre"
ASTRAL_LOG="$LOGS_DIR/astral_run.log"

echo "Running ASTRAL (C++ version via Conda)..."
# We call the native astral command provided by the 'aster' conda package
astral \
    -i "$COMBINED_TREES" \
    -o "$SPECIES_TREE" \
    2>&1 | tee "$ASTRAL_LOG"

echo ""
echo "========================================="
echo "ASTRAL Species Tree Estimation Complete!"
echo "Species Tree saved to: $SPECIES_TREE"
echo "Time: $(date)"
echo "========================================="
