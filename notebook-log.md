## Project Title: Mannose-associated conserved genes in 267 yeast species (CDS/proteins → MSA)
**Notebook date:** 2026-02-10
### Goal:
### Use orthogroup-defined gene families (from Orthofinder across 1,154 Y1000+ species) and restrict them to your 267 focal species, then align each orthogroup separately and trim low-quality regions to enable downstream conservation analysis and phylogenetic inference. The main aim is to analyze **267 yeast species** from diverse taxa to identify **genes conserved in mannose-utilizing species** compared with species that do not metabolize mannose.

### Core logic:(CDS/proteins → MSA → trimming → gene trees)

1. My ML models identify orthogroups (OGs) associated with mannose phenotype.
2. An OG is a gene-family “bucket”, it may contain a sinfle gene though. 
3. Each species can contribute 0, 1, or many genes to a given OG (paralogs happen).Therefore for this analyses, MSA is done per orthogroup, not “per species”.
Each species can contribute 0, 1, or many genes to an OG (gene family / paralogs).

### Core workflow (per orthogroup)
1. QC input sequences (CDS/proteins/genomes as available)
2. Build comparable gene sets across species (Orthogroups / for curated genes)
3. **Subset** full Y1000+ orthogroup FASTAs (1,154 species) down to **my 267 species**
4. Run **MSA per orthogroup** (using MAFFT with the following configurations: for OGs that were small ie <200 sequences, and also skipped OGs with <3 sequecnes. For large OGs,is chose  --auto selected:FFT-NS-2 (fast progressive method, guide trees built twice) but for small OGs with < with 3 sequences,  used L-INS-i (--localpair --maxiterate 1000) = highly accurate ). I did this because accuracy is maximized where feasible (small OGs) and computation remains feasible for massive gene families (large OGs)
5. Trim alignments (TrimAl)
6. Infer gene trees / model selection / support (IQ-TREE)
7. Downstream: quantify conservation / enrichment between phenotype groups


### Environment setup (conda)

I created a working environment to preevent tools interferring with other system tools:

```bash
conda create -n phylo_env -c bioconda -c conda-forge iqtree mafft trimal fasttree
conda activate phylo_env
```

```python
# (In the terminal, I ran):
!conda --version
# Optional: verify tools are visible
!which mafft || true
!which trimal || true
!which iqtree3 || true
!which iqtree || true
```

## Assumed repository layout

I ran all analyses from the root repository root ( `Phylogenetics-repository-for-PLANTPATH-563`).  
All paths below are relative to that root.

### Inputs I already have

- `data/processed/species_267.txt` : list of the 267 species IDs (one per line)
- `data/processed/03_ml_orthogroups/01_og_lists/orthogroup_list.txt` : **master OG list** (union of my 6 ML models which contained 9 OGs)
- Orthofinder OG FASTAs (all species):  
  `data/raw/y1000p_orthofinder/Orthogroup_Sequences/OGxxxxxxx.fasta`

Already created:
- Subset OG FASTAs to your 267 species:  
  `data/processed/03_ml_orthogroups/05_og_fastas_267/*.fasta`
- MAFFT alignments per OG:  
  `data/processed/03_ml_orthogroups/06_alignments/*.aln.fasta`

### Created and verified directories

```python
!mkdir -p data/processed/03_ml_orthogroups/logs
!mkdir -p data/processed/03_ml_orthogroups/05_og_fastas_267
!mkdir -p data/processed/03_ml_orthogroups/06_alignments
!mkdir -p data/processed/03_ml_orthogroups/07_trimmed_alignments
!mkdir -p data/processed/03_ml_orthogroups/08_gene_trees
!ls -lah data/processed/03_ml_orthogroups | sed -n '1,120p'
```

### Step A — Subset Orthogroup FASTAs to 267 species (seqkit)

This step reads each OG in my `orthogroup_list.txt`, opens the full OG FASTA from Orthofinder (1,154 species),
and filters sequences whose headers contain **any** species ID from `species_267.txt`.

**Important header pattern note:** My OG FASTA headers look like (this initally creates parsing problems but I figured it out and solved it):
`>g009075.m1|aspergillus_nidulans.final`  
So my species list should contain the `*.final` IDs (matching the portion after the pipe).

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

```python
# Run
!./og_fastas_267.sh

# Count outputs
!ls -1 data/processed/03_ml_orthogroups/05_og_fastas_267/*.fasta 2>/dev/null | wc -l || true
!head -n 20 data/processed/03_ml_orthogroups/logs/og_counts_267.tsv || true
!tail -n 20 data/processed/03_ml_orthogroups/logs/zero_after_filter.log || true
```

### Step B — MAFFT alignments per OG (robust loop)

Safeguards:
- skip OGs with <3 sequences
- use more accurate method for ≤200 sequences
- use `--auto` for large OGs (scalable)
- write a run log

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

### Run MAFFT + verify outputs with the following parameters:

```python
!mafft --version | head -n 2
!./mafft_run.sh
!ls -1 data/processed/03_ml_orthogroups/06_alignments/*.aln.fasta 2>/dev/null | wc -l || true
!tail -n 30 data/processed/03_ml_orthogroups/logs/mafft_run.log || true
```
As already indicated, I set FFT-NS-2 with the flag --auto to build trees twice for large OGs and L-INS-i --localpair --maxiterate 1000 for smaller OGs.
### Step C — Trim alignments with TrimAl (robust + avoids `HISTTIMEFORMAT` unbound errors: I ran into this error and it took some time to troubleshoot it)

You hit an interactive shell error:
`HISTTIMEFORMAT: unbound variable`

This happens with `set -u` when something references an unset env var.
The script below defensively sets `HISTTIMEFORMAT` if missing.
It also writes to a `.tmp` file and then moves it into place (safer).

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

  # macOS: time -l prints memory/time stats; if unavailable, just run trimal
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

I trimmed each MAFFT alignment with trimal -automated1 because large OGs were aligned with FFT-NS-2, which is fast but rough. This introduces introduces noisy columns, produces gap-rich regions and 
can inflate false phylogenetic signal. 
Because I was running a batch pipeline across many OGs, with widely varying properties (length, divergence, gap patterns), i chose -automated1 is ideal because it adapts trimming automatically per alignment
and avoids subjective hand-tuning of thresholds. 

```python
!trimal -h | head -n 2
!./trim_align.sh
!ls -1 data/processed/03_ml_orthogroups/07_trimmed_alignments/*.trim.aln.fasta 2>/dev/null | wc -l || true
!tail -n 40 data/processed/03_ml_orthogroups/logs/trimal_run.log || true
```


### Step D — Gene trees with IQ-TREE (IQ-TREE 3)

In my environment, `iqtree2` was not found, but `iqtree3` exists, installed `iqtree3`.
Use **Ultrafast bootstrap** with **1000 replicates** via `-B 1000`.

### Run per trimmed alignment
- `-m MFP` : ModelFinder + tree inference
- `-B 1000` : ultrafast bootstrap
- `-T AUTO` : auto threads
- `-pre` : output prefix per OG

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

```python
# NOTE: this can be compute-heavy depending on alignment sizes so I plan to run  this part on HTC cluster
!./iqtree_run.sh
!ls -lah data/processed/03_ml_orthogroups/08_gene_trees | sed -n '1,120p'
!tail -n 40 data/processed/03_ml_orthogroups/logs/iqtree_run.log || true
```

### Next analysis steps (after gene trees)

At this point I have:
- per-OG trimmed alignments (MSA)
- per-OG gene trees + bootstrap support

### Next steps  (things I am thinking of doing)

1. **Define phenotype groups**
   - `mannose-growers` vs `non-growers` (a table mapping species → phenotype)

2. **Quantify conservation per OG**
   - residue conservation scores per column (Shannon entropy, Jensen-Shannon, etc.)
   - per-OG summary: mean entropy, fraction conserved sites, etc.
   - compare distributions between phenotype groups

3. **Compare sequence features by group**
   - group-specific motifs, indels, domain differences
   - test enrichment of conserved substitutions in growers vs non-growers

4. **Handle paralogs (big OGs)**
   - your largest OGs suggest multi-copy gene families
   - options: keep all copies (gene-family view), or reduce to 1 per species (e.g., longest, best-scoring) for clearer phylogenies

5. **Optional: concatenation / species tree**
   - choose a set of single-copy OGs, concatenate alignments, infer species tree
   - or use gene-tree methods (ASTRAL) if desired

### Reproducibility checklist

- Keep all key lists under version control:
  - `species_267.txt`
  - `orthogroup_list.txt` (master OG list)
- Log every run:
  - `logs/mafft_run.log`
  - `logs/trimal_run.log`
  - `logs/iqtree_run.log`
- Capture tool versions:
  - `mafft --version`
  - `trimal -h`
  - `iqtree3 -v`
