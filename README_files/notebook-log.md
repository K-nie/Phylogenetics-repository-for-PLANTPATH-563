# Lab Notebook — Mannose-associated conserved genes in 267 yeast species

**Notebook date:** 2026-02-10

## Goal

To use orthogroup-defined gene families (from OrthoFinder across 1,154 Y1000+ species) and restrict them to 267 focal species, then align each orthogroup separately and trim low-quality regions to enable downstream conservation analysis and phylogenetic inference. The main aim is to analyze **267 yeast species** from diverse taxa to identify **genes conserved in mannose-utilizing species** compared with species that do not metabolize mannose.

---

## Core Logic (CDS/proteins → MSA → trimming → gene trees)

1. ML models identify orthogroups (OGs) associated with the mannose phenotype.
2. An OG is a gene-family "bucket" — it may contain a single gene, though paralogs are common.
3. Each species can contribute 0, 1, or many genes to a given OG. Therefore, MSA is done per orthogroup, not per species.

---

## Core Workflow (per orthogroup)

1. QC input sequences (CDS/proteins/genomes as available)
2. Build comparable gene sets across species (orthogroups, or curated genes)
3. **Subset** full Y1000+ orthogroup FASTAs (1,154 species) down to **267 focal species**
4. Run **MSA per orthogroup** using MAFFT:
   - OGs with ≤ 200 sequences: `--localpair --maxiterate 1000` (L-INS-i — high accuracy)
   - OGs with > 200 sequences: `--auto` (FFT-NS-2 — fast progressive, guide trees built twice)
   - OGs with < 3 sequences: skipped
   - Rationale: accuracy is maximized where feasible (small OGs); computation remains tractable for large gene families
5. Trim alignments (TrimAl)
6. Infer gene trees, model selection, and bootstrap support (IQ-TREE)
7. Downstream: quantify conservation and enrichment between phenotype groups

---

## Environment Setup (conda)

A dedicated environment was created to prevent tools from interfering with other system tools:

```bash
conda create -n phylo_env -c bioconda -c conda-forge iqtree mafft trimal fasttree
conda activate phylo_env
```

Verify tools are accessible:

```bash
conda --version
which mafft
which trimal
which iqtree3 || which iqtree
```

---

## Repository Layout

All analyses were run from the repository root (`Phylogenetics-repository-for-PLANTPATH-563`). All paths below are relative to that root.

### Inputs

- `data/processed/species_267.txt` — list of 267 species IDs (one per line)
- `data/processed/03_ml_orthogroups/01_og_lists/orthogroup_list.txt` — master OG list (union of 6 ML models, containing 9 OGs)
- OrthoFinder OG FASTAs (all species): `data/raw/y1000p_orthofinder/Orthogroup_Sequences/OGxxxxxxx.fasta`

### Outputs (created during preprocessing)

- Subset OG FASTAs (267 species): `data/processed/03_ml_orthogroups/05_og_fastas_267/*.fasta`
- MAFFT alignments per OG: `data/processed/03_ml_orthogroups/06_alignments/*.aln.fasta`
- TrimAl trimmed alignments: `data/processed/03_ml_orthogroups/07_trimmed_alignments/*.trim.aln.fasta`
- IQ-TREE gene trees (local test): `data/processed/03_ml_orthogroups/08_gene_trees/`

### Create and verify directories

```bash
mkdir -p data/processed/03_ml_orthogroups/logs
mkdir -p data/processed/03_ml_orthogroups/05_og_fastas_267
mkdir -p data/processed/03_ml_orthogroups/06_alignments
mkdir -p data/processed/03_ml_orthogroups/07_trimmed_alignments
mkdir -p data/processed/03_ml_orthogroups/08_gene_trees
ls -lah data/processed/03_ml_orthogroups
```

---

## Step A — Subset Orthogroup FASTAs to 267 Species (seqkit)

This step reads each OG in `orthogroup_list.txt`, opens the full OG FASTA from OrthoFinder (1,154 species), and filters sequences whose headers match any species ID from `species_267.txt`.

**Header pattern note:** OG FASTA headers follow the pattern `>g009075.m1|aspergillus_nidulans.final`. The species list must therefore contain the `*.final` IDs — the portion after the pipe — to match correctly.

### Script: `og_fastas_267.sh`

```bash
cat > og_fastas_267.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

OGDIR="data/raw/y1000p_orthofinder/Orthogroup_Sequences"
OUTDIR="data/processed/03_ml_orthogroups/05_og_fastas_267"
ALOG="data/processed/03_ml_orthogroups/logs"
LIST="data/processed/03_ml_orthogroups/01_og_lists/orthogroup_list.txt"
SP="data/processed/species_267.txt"

mkdir -p "$OUTDIR" "$ALOG"

MISSING_LOG="${ALOG}/missing_og_fastas.log"
COUNTS_TSV="${ALOG}/og_counts_267.tsv"
ZERO_LOG="${ALOG}/zero_after_filter.log"

: > "$MISSING_LOG"
: > "$COUNTS_TSV"
: > "$ZERO_LOG"

command -v seqkit >/dev/null 2>&1 || { echo "[ERROR] seqkit not found in PATH" >&2; exit 1; }

echo -e "OG\tNseq_267" >> "$COUNTS_TSV"

while read -r og; do
  [[ -z "${og// }" ]] && continue
  [[ "$og" =~ ^# ]] && continue

  in="${OGDIR}/${og}.fasta"
  out="${OUTDIR}/${og}.fasta"

  if [ ! -s "$in" ]; then
    echo "[WARN] Missing: $in" | tee -a "$MISSING_LOG" >&2
    continue
  fi

  seqkit grep -f "$SP" "$in" > "$out" || true

  n=$(grep -c '^>' "$out" 2>/dev/null || true)
  echo -e "${og}\t${n}" >> "$COUNTS_TSV"

  if [ "${n}" -eq 0 ]; then
    echo "[INFO] ${og} has 0 sequences after filtering: $out" | tee -a "$ZERO_LOG" >&2
  fi
done < "$LIST"

echo "[DONE] Subsetting complete."
echo "  Output FASTAs : $OUTDIR"
echo "  Counts TSV    : $COUNTS_TSV"
echo "  Missing OGs   : $MISSING_LOG"
echo "  Zero after    : $ZERO_LOG"
EOF

chmod +x og_fastas_267.sh
```

### Run subsetting + sanity checks

```bash
./og_fastas_267.sh

# Count outputs
ls -1 data/processed/03_ml_orthogroups/05_og_fastas_267/*.fasta 2>/dev/null | wc -l
head -n 20 data/processed/03_ml_orthogroups/logs/og_counts_267.tsv
tail -n 20 data/processed/03_ml_orthogroups/logs/zero_after_filter.log
```

---

## Step B — Multiple Sequence Alignment (MAFFT)

Safeguards built into the script:
- Skip OGs with fewer than 3 sequences
- Use L-INS-i (`--localpair --maxiterate 1000`) for OGs with ≤ 200 sequences
- Use `--auto` (FFT-NS-2) for OGs with > 200 sequences
- Write a run log for reproducibility

### Script: `mafft_run.sh`

```bash
cat > mafft_run.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

IN="data/processed/03_ml_orthogroups/05_og_fastas_267"
OUT="data/processed/03_ml_orthogroups/06_alignments"
LOG="data/processed/03_ml_orthogroups/logs/mafft_run.log"

mkdir -p "$OUT" "$(dirname "$LOG")"
: > "$LOG"

command -v mafft >/dev/null 2>&1 || { echo "[ERROR] mafft not found in PATH" >&2; exit 1; }

shopt -s nullglob
files=("$IN"/*.fasta)
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
  echo "[ERROR] No OG FASTAs found: $IN/*.fasta" | tee -a "$LOG" >&2
  exit 1
fi

for f in "${files[@]}"; do
  og=$(basename "$f" .fasta)
  n=$(grep -c '^>' "$f" || true)

  if [ "$n" -lt 3 ]; then
    echo "[SKIP] $og : only $n sequences" | tee -a "$LOG"
    continue
  fi

  if [ "$n" -le 200 ]; then
    algo="--maxiterate 1000 --localpair"
  else
    algo="--auto"
  fi

  echo "[RUN] $og : $n sequences : mafft $algo" | tee -a "$LOG"
  mafft $algo "$f" > "${OUT}/${og}.aln.fasta"
done

echo "[DONE] MAFFT complete: wrote alignments to $OUT" | tee -a "$LOG"
EOF

chmod +x mafft_run.sh
```

### Run MAFFT + verify outputs

```bash
mafft --version | head -n 2
./mafft_run.sh
ls -1 data/processed/03_ml_orthogroups/06_alignments/*.aln.fasta 2>/dev/null | wc -l
tail -n 30 data/processed/03_ml_orthogroups/logs/mafft_run.log
```

**MAFFT mode summary:** `--auto` (FFT-NS-2) was used for large OGs — this builds guide trees twice for improved accuracy at scale. L-INS-i (`--localpair --maxiterate 1000`) was used for smaller OGs where accuracy is maximized without excessive computation.

---

## Step C — Alignment Trimming (TrimAl)

Each MAFFT alignment was trimmed with `trimal -automated1`.

**Rationale:** Large OGs aligned with FFT-NS-2 can introduce noisy, gap-rich columns that inflate false phylogenetic signal. The `-automated1` flag adapts trimming thresholds automatically per alignment, making it ideal for a batch pipeline across OGs with widely varying length, divergence, and gap patterns — no manual threshold tuning required.

**Troubleshooting note:** An `HISTTIMEFORMAT: unbound variable` error was encountered when running with `set -u`. This occurs because interactive shells or plugins reference this variable without setting it. The script defensively initializes it with `: "${HISTTIMEFORMAT:=}"`. Output is also written to a `.tmp` file before being moved into place to prevent partial outputs from corrupting the pipeline.

### Script: `trim_align.sh`

```bash
cat > trim_align.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

# Avoid "unbound variable" crashes from shells/plugins referencing HISTTIMEFORMAT
: "${HISTTIMEFORMAT:=}"

IN="data/processed/03_ml_orthogroups/06_alignments"
OUT="data/processed/03_ml_orthogroups/07_trimmed_alignments"
LOG="data/processed/03_ml_orthogroups/logs/trimal_run.log"

mkdir -p "$OUT" "$(dirname "$LOG")"
: > "$LOG"

command -v trimal >/dev/null 2>&1 || { echo "[ERROR] trimal not found in PATH" | tee -a "$LOG" >&2; exit 1; }

shopt -s nullglob
files=("$IN"/*.aln.fasta)
shopt -u nullglob

echo "[INFO] Inputs matched: ${#files[@]} files" | tee -a "$LOG"
if [ ${#files[@]} -eq 0 ]; then
  echo "[ERROR] No MAFFT alignments found: $IN/*.aln.fasta" | tee -a "$LOG" >&2
  exit 1
fi

for aln in "${files[@]}"; do
  og=$(basename "$aln" .aln.fasta)
  out="${OUT}/${og}.trim.aln.fasta"
  tmp="${out}.tmp"

  echo "[TRIMAL] $og" | tee -a "$LOG"

  if command -v /usr/bin/time >/dev/null 2>&1; then
    /usr/bin/time -l trimal -in "$aln" -out "$tmp" -automated1 >>"$LOG" 2>&1
  else
    trimal -in "$aln" -out "$tmp" -automated1 >>"$LOG" 2>&1
  fi

  mv -f "$tmp" "$out"

  if [ ! -s "$out" ] || ! grep -q '^>' "$out"; then
    echo "[ERROR] Output empty for $og: $out" | tee -a "$LOG" >&2
    exit 1
  fi
done

echo "[DONE] TrimAl complete: wrote ${#files[@]} files to $OUT" | tee -a "$LOG"
EOF

chmod +x trim_align.sh
```

### Run TrimAl + verify

```bash
trimal -h | head -n 2
./trim_align.sh
ls -1 data/processed/03_ml_orthogroups/07_trimmed_alignments/*.trim.aln.fasta 2>/dev/null | wc -l
tail -n 40 data/processed/03_ml_orthogroups/logs/trimal_run.log
```

---

## Step D — Gene Trees with IQ-TREE (local test)

IQ-TREE 3 (`iqtree3`) is used. The script falls back to `iqtree` if `iqtree3` is not found. Ultrafast bootstrap is run with 1,000 replicates (`-B 1000`).

**Note:** This local run is used to validate the preprocessing pipeline. The full model selection and batch tree inference for all 9 OGs runs on the GLBRC HT Condor cluster (see main README for cluster pipeline details).

### Script: `iqtree_run.sh`

```bash
cat > iqtree_run.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

IN="data/processed/03_ml_orthogroups/07_trimmed_alignments"
OUT="data/processed/03_ml_orthogroups/08_gene_trees"
LOG="data/processed/03_ml_orthogroups/logs/iqtree_run.log"

mkdir -p "$OUT" "$(dirname "$LOG")"
: > "$LOG"

# Prefer iqtree3, fallback to iqtree
IQ=""
if command -v iqtree3 >/dev/null 2>&1; then
  IQ="iqtree3"
elif command -v iqtree >/dev/null 2>&1; then
  IQ="iqtree"
else
  echo "[ERROR] iqtree3/iqtree not found in PATH" | tee -a "$LOG" >&2
  exit 1
fi

shopt -s nullglob
files=("$IN"/*.trim.aln.fasta)
shopt -u nullglob

echo "[INFO] Inputs matched: ${#files[@]} files" | tee -a "$LOG"
if [ ${#files[@]} -eq 0 ]; then
  echo "[ERROR] No trimmed alignments found: $IN/*.trim.aln.fasta" | tee -a "$LOG" >&2
  exit 1
fi

for aln in "${files[@]}"; do
  og=$(basename "$aln" .trim.aln.fasta)
  pre="${OUT}/${og}"

  echo "[RUN] $og : $IQ -m MFP -B 1000" | tee -a "$LOG"
  $IQ -s "$aln" -m MFP -B 1000 -T AUTO -pre "$pre" >>"$LOG" 2>&1
done

echo "[DONE] IQ-TREE complete: wrote outputs to $OUT" | tee -a "$LOG"
EOF

chmod +x iqtree_run.sh
```

### Run IQ-TREE + verify outputs

```bash
# Note: compute-heavy for large alignments — run on HTC cluster for full batch
./iqtree_run.sh
ls -lah data/processed/03_ml_orthogroups/08_gene_trees
tail -n 40 data/processed/03_ml_orthogroups/logs/iqtree_run.log
```

**IQ-TREE parameters:**

| Flag | Purpose |
|------|---------|
| `-m MFP` | ModelFinder + tree inference in one run |
| `-B 1000` | Ultrafast bootstrap (1,000 replicates) |
| `-T AUTO` | Automatic thread detection |
| `-pre` | Output prefix per OG |

---

## Next Analysis Steps (after gene trees)

At this point the pipeline has produced:
- Per-OG trimmed alignments (MSA)
- Per-OG gene trees with bootstrap support

Planned downstream analyses:

1. **Define phenotype groups** — map species to mannose-grower vs. non-grower phenotype
2. **Quantify conservation per OG** — residue conservation scores per column (Shannon entropy, Jensen-Shannon divergence); compare distributions between phenotype groups
3. **Compare sequence features by group** — group-specific motifs, indels, domain differences; enrichment of conserved substitutions in growers vs. non-growers
4. **Handle paralogs in large OGs** — options: keep all copies (gene-family view), or reduce to one per species (e.g., longest/best-scoring) for cleaner phylogenies
5. **Species tree (optional)** — concatenate single-copy OGs and infer a species tree, or use gene-tree coalescence methods (ASTRAL)

---

## Reproducibility Checklist

Keep all key input lists under version control:
- `species_267.txt`
- `orthogroup_list.txt` (master OG list)

Log every run:
- `logs/mafft_run.log`
- `logs/trimal_run.log`
- `logs/iqtree_run.log`

Capture and record tool versions:
```bash
mafft --version
trimal -h | head -n 2
iqtree3 -v
```

---

## Genome Alignment — pygenomeviz / progressiveMauve

### Chosen Method and Rationale

**Method:** pygenomeviz (wrapping progressiveMauve from the Mauve genome alignment package, Darling et al., 2010).

**Rationale:** In the Alignathon benchmark, progressiveMauve produced high-quality, transitively closed whole-genome alignments (consistent pairwise homologies without reference bias) and performed well on simulated datasets with rearrangements. It was highlighted for excellent precision in certain contexts (e.g., primate-like data) and for handling gene gain/loss and structural variation effectively. It suits this project's dataset — species from many clades with expected significant sequence differences, inversions, rearrangements, and variable gene content. progressiveMauve's locally collinear block (LCB) approach and progressive guide-tree strategy make it effective for detecting conserved regions across subsets of taxa while tolerating structural differences.

**Scope:** For the homework, only pygenomeviz/progressiveMauve is being run. For the final project, at least two methods will be compared (e.g., pygenomeviz + Cactus) to evaluate alignment quality across metrics such as coverage of conserved regions, handling of rearrangements, and downstream effects on phylogeny or variant detection.

### Algorithm Description

pygenomeviz/progressiveMauve is a progressive multiple genome aligner. It identifies maximal unique matches (anchors/MUMs) to define locally collinear blocks (LCBs — regions of conserved sequence order), builds a phylogenetic guide tree from shared gene content, aligns sequences progressively along the tree, refines anchors recursively, performs gapped alignment within LCBs, and uses a homology hidden Markov model (HMM) to filter spurious matches. It also optimizes a sum-of-pairs breakpoint score to accurately detect rearrangement breakpoints even with unequal gene content.

### Assumptions

- Genomes share orthologous regions detectable via unique anchors
- Large collinear blocks exist despite rearrangements and indels
- Divergence level allows reliable anchor detection (typically effective at > ~50–70% identity in conserved regions)
- A reasonable phylogenetic structure exists for the guide tree

### Limitations

- Computationally demanding for large numbers of genomes or very large eukaryotic genomes (runtime scales roughly with the cube of sequence number in worst cases)
- May fragment LCBs or miss alignments in extremely divergent regions (< 50% identity) or regions with heavy repeats/duplications
- High memory usage is possible
- Breakpoint detection may require parameter tuning for highly rearranged data
- Less scalable than modern graph-based tools (e.g., Cactus) for massive or highly complex variation — pygenomeviz's updated implementation partially addresses this
