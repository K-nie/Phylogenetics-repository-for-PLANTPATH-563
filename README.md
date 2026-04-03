# Phylogenetics Repository — PLANTPATH 563

**Mannose-associated conserved genes in 267 yeast species**

> Phylogenomic pipeline to identify and analyze gene families associated with mannose utilization across 267 *Saccharomycotina* yeast species using orthogroups, multiple sequence alignment, trimming, and gene tree inference.

**Author:** Benjamin Narh-Madey | **Cluster:** GLBRC scarcity-ap-1 | **Conda env:** `phylo_env` | **Last updated:** March 2026

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Research Questions](#2-research-questions)
3. [Dataset](#3-dataset)
4. [Repository Structure](#4-repository-structure)
5. [Environment Setup](#5-environment-setup)
6. [Preprocessing Pipeline (Local)](#6-preprocessing-pipeline-local)
7. [Analysis Pipeline (Cluster)](#7-analysis-pipeline-cluster)
8. [Scripts Reference](#8-scripts-reference)
9. [Output Files](#9-output-files)
10. [Tools & Environment](#10-tools--environment)
11. [Data Formats & File Specifications](#11-data-formats--file-specifications)
12. [Functions Reference](#12-functions-reference)
13. [D-statistic Interpretation](#13-d-statistic-interpretation)
14. [Known Issues & Resolutions](#14-known-issues--resolutions)
15. [Current Status & Next Steps](#15-current-status--next-steps)

---

## 1. Project Overview

This project analyzes 267 focal yeast species drawn from the Y1000+ dataset to identify genes conserved in mannose-utilizing species compared with non-utilizers. Machine learning models identified orthogroups (OGs) associated with the mannose growth phenotype. Each orthogroup represents a gene family; each species may contribute 0, 1, or multiple genes per OG (paralogs are expected). Therefore, MSA is performed per orthogroup — not per species.

The full pipeline spans two stages:

- **Preprocessing (local):** Subset the full Y1000+ OrthoFinder FASTAs (1,154 species) down to 267 focal species, align each orthogroup, trim low-quality columns, and run a local IQ-TREE test.
- **Cluster analysis (GLBRC HT Condor):** Parallel ModelFinder model selection, maximum likelihood trees (RAxML-NG), Bayesian trees (BEAST2), and D-statistic phylogenetic signal analysis (R).

A test run was completed on OG0000007 to validate the full pipeline before batch submission. Best-fit model: **Q.YEAST+F+R10 (BIC)**.

---

## 2. Research Questions

1. **How has the genetic basis for mannose metabolism evolved over the *Saccharomycotina* subphylum?**
   By generating Maximum Likelihood trees (RAxML-NG) and time-calibrated Bayesian trees (BEAST2), we build an evolutionary roadmap. MRCA anchors allow us to see not just how mannose metabolism genes diverged, but approximately *when* (mya) those splits occurred across 267 species.

2. **Does the evolutionary history of these genes dictate a species' ability to grow on mannose?**
   By mapping binary trait data (`1`=growth, `0`=no growth) from `mannose_traits.csv` onto BEAST2 chronograms using FigTree, we can pinpoint the evolutionary eras where this metabolic ability was gained, lost, or conserved.

3. **Does mannose growth ability show phylogenetic signal across yeast gene trees?**
   The D-statistic script (`run_dstatistic.R`) calculates whether the mannose growth trait is strongly linked to species evolutionary history or randomly distributed across the subphylum.

4. **Which substitution model best fits each orthogroup alignment?**
   IQ-TREE ModelFinder tests 1,232 models per alignment and selects the best by BIC. The selected model is passed directly to RAxML-NG, replacing the old pipeline that hardcoded `LG` for all orthogroups.

---

## 3. Dataset

| Item | Description |
|------|-------------|
| Orthogroup alignments | 9 trimmed protein alignment files (`OG*.trim.aln.fasta`) in `data/` |
| Species count | 267 yeast species from the Y1000+ *Saccharomycotina* dataset |
| Trait data | Binary mannose growth data (`1`=growth, `0`=no growth) in `data/mannose_traits.csv` |
| Source archives | `00_downloads/y1000p_orthofinder.tar.gz` and `y1000p_pep_files.tar.gz` |
| Raw OG FASTAs | `data/raw/y1000p_orthofinder/Orthogroup_Sequences/` (1,154 species, all OGs) |

> **Note on OG0007179:** This orthogroup was removed from the Condor analysis. The outgroup (*Lipomyces starkeyi*) is absent, it contains fewer than 4 genes, and it cannot be properly rooted for D-statistic calculation.

---

## 4. Repository Structure

```
Phylogenetics-repository-for-PLANTPATH-563/
├── data/
│   ├── OG*.trim.aln.fasta            # Input alignments for cluster pipeline
│   ├── mannose_traits.csv            # Species × mannose growth trait (267 rows)
│   ├── raw/
│   │   └── y1000p_orthofinder/
│   │       └── Orthogroup_Sequences/ # Full OG FASTAs (1,154 species)
│   └── processed/
│       ├── species_267.txt           # List of 267 focal species
│       └── 03_ml_orthogroups/
│           ├── 01_og_lists/
│           │   └── orthogroup_list.txt
│           ├── 05_og_fastas_267/     # OG FASTAs subset to 267 species
│           ├── 06_alignments/        # MAFFT outputs
│           ├── 07_trimmed_alignments/ # TrimAl outputs
│           ├── 08_gene_trees/        # IQ-TREE outputs (local test)
│           └── logs/                 # Preprocessing logs
├── beast_xml/                        # BEAST2 XML configs (one per OG)
├── iqtree_logs/                      # IQ-TREE ModelFinder output per OG
├── results/
│   ├── OG*.raxml.*                   # RAxML-NG output trees
│   └── beast/                        # BEAST2 output trees
├── logs/                             # Condor + runtime logs
├── scripts/
│   ├── modelfinder_jobs.sub          # Phase 1 — Condor submit file
│   ├── run_modelfinder_job.sh        # Phase 1 — per-job ModelFinder script
│   ├── collect_models.sh             # Phase 1 — collects models post-run
│   ├── run_modelfinder.sh            # Phase 1 — legacy serial version (backup)
│   ├── raxml_jobs.sub                # Phase 2a — Condor submit file
│   ├── run_raxml_job.sh              # Phase 2a — per-job RAxML-NG script
│   ├── beast_jobs.sub                # Phase 2b — Condor submit file
│   ├── run_beast_job.sh              # Phase 2b — per-job BEAST2 script
│   ├── generate_beast_xml.py         # Phase 2b — generates XML configs
│   ├── run_dstatistic.R              # Phase 3 — D-statistic analysis
│   ├── modelfinder_jobs.txt          # Generated: ls data/OG* | xargs basename
│   ├── jobs_with_models.txt          # Generated: by collect_models.sh
│   └── beast_jobs.txt                # Generated: ls beast_xml/*.xml | ...
├── bin/                              # raxml-ng binary
├── 00_downloads/                     # Raw downloaded data
├── README_files/                     # Source README documents (archived)
└── phylo_env.tar.gz                  # Packaged conda environment (for Condor nodes)
```

---

## 5. Environment Setup

All tools are managed in a single conda environment named `phylo_env`.

```bash
conda create -n phylo_env -c bioconda -c conda-forge \
    iqtree mafft trimal fasttree seqkit raxml-ng beast2 beagle-lib \
    r-base r-caper r-picante r-phytools r-ape

conda activate phylo_env
```

To verify tools are accessible:

```bash
conda activate phylo_env
which mafft && mafft --version
which trimal && trimal -h | head -n 2
which iqtree3 || which iqtree
which raxml-ng && raxml-ng --version
```

To repackage the environment for deployment to Condor execute nodes:

```bash
conda pack -n phylo_env -o phylo_env.tar.gz
```

---

## 6. Preprocessing Pipeline (Local)

This stage runs locally (or on the submit node) before submitting cluster jobs. It produces the trimmed orthogroup alignments (`OG*.trim.aln.fasta`) used by all downstream analyses.

### Core Logic

- ML models identify OGs associated with mannose phenotype
- An OG is a gene-family "bucket" — each species contributes 0, 1, or many genes (paralogs)
- MSA is done **per orthogroup**, not per species
- Preprocessing flow: subset → align (MAFFT) → trim (TrimAl) → local gene tree test (IQ-TREE)

### Step A — Subset Orthogroup FASTAs to 267 Species

Each OG FASTA from OrthoFinder (1,154 species) is filtered to retain only sequences from the 267 focal species using `seqkit grep`.

**Important:** OG FASTA headers follow the pattern `>g009075.m1|aspergillus_nidulans.final`, so `species_267.txt` must contain the `*.final` IDs matching the portion after the pipe.

Script: `og_fastas_267.sh`

```bash
# Run subsetting
./og_fastas_267.sh

# Verify outputs
ls -1 data/processed/03_ml_orthogroups/05_og_fastas_267/*.fasta | wc -l
head -n 20 data/processed/03_ml_orthogroups/logs/og_counts_267.tsv
```

### Step B — Multiple Sequence Alignment (MAFFT)

Each orthogroup is aligned independently. OGs with fewer than 3 sequences are skipped.

| OG size | MAFFT mode | Reason |
|---------|------------|--------|
| ≤ 200 sequences | `--localpair --maxiterate 1000` (L-INS-i) | High accuracy |
| > 200 sequences | `--auto` (FFT-NS-2) | Scalable; builds guide tree twice |

Script: `mafft_run.sh`

```bash
./mafft_run.sh
ls -1 data/processed/03_ml_orthogroups/06_alignments/*.aln.fasta | wc -l
tail -n 30 data/processed/03_ml_orthogroups/logs/mafft_run.log
```

### Step C — Alignment Trimming (TrimAl)

Poorly aligned columns are removed using `trimal -automated1`. This method was chosen because:
- Large OGs aligned with FFT-NS-2 can introduce noisy, gap-rich columns that inflate false phylogenetic signal
- `-automated1` adapts trimming thresholds automatically per alignment, avoiding manual parameter tuning across diverse OG sizes

Script: `trim_align.sh`

```bash
./trim_align.sh
ls -1 data/processed/03_ml_orthogroups/07_trimmed_alignments/*.trim.aln.fasta | wc -l
tail -n 40 data/processed/03_ml_orthogroups/logs/trimal_run.log
```

> **Troubleshooting note:** If you encounter `HISTTIMEFORMAT: unbound variable` errors, the script defensively sets `HISTTIMEFORMAT` to avoid crashes from interactive shell hooks.

### Step D — Local Gene Tree Test (IQ-TREE)

A local IQ-TREE run is used to validate the preprocessing pipeline. The full model selection and tree search for all 9 OGs runs on the cluster (see Section 7).

```bash
# Local test run
./iqtree_run.sh
ls -lah data/processed/03_ml_orthogroups/08_gene_trees
tail -n 40 data/processed/03_ml_orthogroups/logs/iqtree_run.log
```

Parameters used:

| Parameter | Purpose |
|-----------|---------|
| `-m MFP` | Automatic model selection + tree inference |
| `-B 1000` | Ultrafast bootstrap (1,000 replicates) |
| `-T AUTO` | Automatic threading |
| `-pre` | Output prefix per OG |

---

## 7. Analysis Pipeline (Cluster)

The cluster pipeline runs on the **GLBRC HT Condor cluster** (`scarcity-ap-1`). Phases 2a and 2b run in parallel after Phase 1 completes.

```
Phase 1 — Model Selection  (Condor, 9 parallel jobs, ~1.5 hrs total)
  modelfinder_jobs.sub → run_modelfinder_job.sh → collect_models.sh → jobs_with_models.txt

Phase 2a — ML Trees  (Condor, 9 parallel jobs, ~2 hrs each)
  raxml_jobs.sub → run_raxml_job.sh → results/OG*.raxml.support

Phase 2b — Bayesian Trees  (Condor, 9 parallel jobs, ~24 hrs each)
  beast_jobs.sub → run_beast_job.sh → results/beast/OG*.beast.tree

Phase 3 — D-statistic  (submit node, after Phase 2 complete)
  run_dstatistic.R → results/dstatistic_results.tsv
```

### Phase 1 — Model Selection (Condor, parallel)

Phase 1 runs all 9 orthogroups in parallel on Condor execute nodes (~1.5 hrs vs. ~13 hrs for the legacy serial approach).

```bash
conda activate phylo_env

# Step 1: Generate alignment job list
ls data/OG*.trim.aln.fasta | xargs -n1 basename > scripts/modelfinder_jobs.txt
cat scripts/modelfinder_jobs.txt   # verify all 9 OGs listed

# Step 2: Submit all 9 ModelFinder jobs simultaneously
condor_submit scripts/modelfinder_jobs.sub

# Step 3: Monitor until complete
condor_q

# Step 4: Collect models and generate jobs_with_models.txt
bash scripts/collect_models.sh

# Step 5: Verify output
cat scripts/jobs_with_models.txt
```

When complete, `jobs_with_models.txt` will contain one line per orthogroup (tab-separated):

```
OG0000007.trim.aln.fasta    Q.YEAST+F+R10
OG0000008.trim.aln.fasta    LG+F+R10
...
```

**Completed model selection results:**

| Orthogroup | Best Model (BIC) | Model Components |
|------------|-----------------|-----------------|
| OG0000007 | Q.YEAST+F+R10 | Yeast-specific matrix + empirical freqs + 10 FreeRate categories |
| OG0000008 | LG+F+R10 | General LG matrix + empirical freqs + 10 FreeRate categories |
| OG0000009 | LG+R10 | General LG matrix + 10 FreeRate categories |
| OG0000020 | Q.YEAST+R10 | Yeast-specific matrix + 10 FreeRate categories |
| OG0000029 | LG+F+R9 | General LG matrix + empirical freqs + 9 FreeRate categories |
| OG0000032 | LG+F+R9 | General LG matrix + empirical freqs + 9 FreeRate categories |
| OG0000044 | LG+F+R10 | General LG matrix + empirical freqs + 10 FreeRate categories |
| OG0000067 | LG+F+R9 | General LG matrix + empirical freqs + 9 FreeRate categories |
| OG0007179 | **REMOVED** | No outgroup (*Lipomyces starkeyi* absent), <4 genes, cannot be rooted |

### Phase 2a — Maximum Likelihood Trees (RAxML-NG)

Submit all 8 RAxML jobs simultaneously after Phase 1 completes. Each job runs on 10 CPUs with 8 GB RAM.

```bash
# Verify job list
cat scripts/jobs_with_models.txt

# Submit all ML jobs at once
condor_submit scripts/raxml_jobs.sub

# Monitor
condor_q
```

**RAxML-NG settings:** 10 parsimony + 10 random starting trees, 100 bootstrap replicates, per-OG substitution model from `jobs_with_models.txt`, seed `12345`.

### Phase 2b — Bayesian Inference (BEAST2)

Runs in parallel with Phase 2a. Settings: relaxed lognormal clock, Birth-Death tree prior, 50 million MCMC steps.

```bash
# Step 1: Generate one XML config per orthogroup
python3 scripts/generate_beast_xml.py

# Step 2: Inject MRCA Calibration Priors (manual — see table below)

# Step 3: Generate job list
ls beast_xml/*.xml | xargs -n1 basename | sed 's/\.xml//' > scripts/beast_jobs.txt

# Step 4: Submit all BEAST2 jobs
condor_submit scripts/beast_jobs.sub

# After jobs complete — summarize posterior trees (10% burnin)
for f in results/beast/*.trees; do
  og=$(basename $f .trees)
  treeannotator -burnin 10 $f results/beast/${og}.beast.tree
done
```

**MRCA Calibration Priors** (from Shen et al. 2018, MCMCTree Table S2):

| Anchor | Taxa pair | Mean (Ma) | Bounds (Ma) |
|--------|-----------|-----------|-------------|
| Deep | *Babjeviella inositovora* & *Pachysolen tannophilus* | 224.70 | 190.09 – 265.39 |
| Medium | *Zygosaccharomyces rouxii* & *Z. kombuchaensis* | 36.90 | 22.98 – 51.96 |
| Young | *Kazachstania siamensis* & *K. unispora* | 17.27 | 12.36 – 23.83 |

To inject priors: open each `beast_xml/*.xml` in a text editor, locate the prior distribution block, and paste the XML code for each anchor below it. Update `YOUR_TREE_NAME` to match the OG identifier found in the `Tree.t` tag.

> **⚠ IMPORTANT — Convergence check:** Before using Bayesian trees, load the `.log` files into Tracer. All ESS values **must be > 200** (ideally > 500). If not, increase chain length and rerun. Do not proceed to tree summarization until convergence is confirmed.

**Visualization after BEAST2:**

- **FigTree** (publication-ready chronograms) — Open `${og}.beast.tree`. Enable *Node Labels* (set to `posterior` for support or `height` for Ma ages) and *Node Bars* (set to `height_95%_HPD` for age uncertainty).
- **DensiTree** (topological uncertainty) — Open the raw `.trees` file (NOT the summarized tree). This overlays all sampled trees showing regions of high consensus vs. uncertainty.

### Phase 3 — D-statistic Analysis

Run after Phase 2 trees are available. Calculates phylogenetic signal for the mannose growth trait across gene trees.

```bash
# Place mannose_traits.csv in data/ first, then:
Rscript scripts/run_dstatistic.R \
  --trait data/mannose_traits.csv \
  --treedir results/ \
  --output results/dstatistic_results.tsv
```

`mannose_traits.csv` must have three columns: `Species`, `Species_name_in_OGs` (matching tree tip labels exactly), and `Mannose_Class` (0 or 1).

---

## 8. Scripts Reference

| Script | Phase | Description |
|--------|-------|-------------|
| `modelfinder_jobs.sub` | Phase 1 | Condor submit file — queues 9 ModelFinder jobs in parallel |
| `run_modelfinder_job.sh` | Phase 1 | Per-job script: unpacks env, runs IQ-TREE ModelFinder, copies logs |
| `collect_models.sh` | Phase 1 | Scans `iqtree_logs/` and writes `jobs_with_models.txt` after all jobs finish |
| `run_modelfinder.sh` | Phase 1 (legacy) | Original serial version — replaced by Condor submit for cluster use |
| `raxml_jobs.sub` | Phase 2a | Condor submit file — queues 9 RAxML-NG ML tree search jobs in parallel |
| `run_raxml_job.sh` | Phase 2a | Per-job script: unpacks env, runs RAxML-NG, copies results |
| `beast_jobs.sub` | Phase 2b | Condor submit file — queues 9 BEAST2 Bayesian inference jobs in parallel |
| `run_beast_job.sh` | Phase 2b | Per-job script: unpacks env, runs BEAST2, copies results |
| `generate_beast_xml.py` | Phase 2b | Generates one BEAST2 XML config per orthogroup alignment |
| `run_dstatistic.R` | Phase 3 | Calculates D-statistic for mannose trait on gene trees |
| `og_fastas_267.sh` | Preprocessing | Subsets full OG FASTAs to 267 focal species using seqkit |
| `mafft_run.sh` | Preprocessing | Runs MAFFT per OG (L-INS-i for ≤200, --auto for >200 sequences) |
| `trim_align.sh` | Preprocessing | Trims MAFFT alignments with `trimal -automated1` |
| `iqtree_run.sh` | Preprocessing | Local IQ-TREE test run on trimmed alignments |

---

## 9. Output Files

| Directory | File Pattern | Phase | Description |
|-----------|-------------|-------|-------------|
| `iqtree_logs/` | `OG*.trim.aln_model.log` | Phase 1 | IQ-TREE ModelFinder full log per orthogroup |
| `scripts/` | `jobs_with_models.txt` | Phase 1 | Tab-separated: filename + best model (input to Phase 2) |
| `results/` | `OG*.raxml.bestTree` | Phase 2a | Best ML tree topology (no support values) |
| `results/` | `OG*.raxml.support` | Phase 2a | **MAIN RESULT** — Best ML tree with bootstrap support values |
| `results/` | `OG*.raxml.bestModel` | Phase 2a | Best-fit substitution model used by RAxML |
| `results/` | `OG*.raxml.bootstraps` | Phase 2a | All 100 individual bootstrap trees |
| `results/beast/` | `OG*.trace.log` | Phase 2b | BEAST2 MCMC trace — check ESS > 200 in Tracer |
| `results/beast/` | `OG*.trees` | Phase 2b | Posterior tree distribution from BEAST2 |
| `results/beast/` | `OG*.beast.tree` | Phase 2b | **MAIN RESULT** — Summarized Bayesian tree (after TreeAnnotator) |
| `results/` | `dstatistic_results.tsv` | Phase 3 | D-statistic per orthogroup: D value, p-values, interpretation |
| `logs/` | `*.out / *.err` | All | Condor stdout/stderr logs per job |

---

## 10. Tools & Environment

All tools are installed in the `phylo_env` conda environment at `/home/glbrc.org/narhmadey/.conda/envs/phylo_env`. The environment is also packaged as `phylo_env.tar.gz` for deployment to Condor execute nodes.

| Tool | Version | Channel | Purpose | Status |
|------|---------|---------|---------|--------|
| IQ-TREE | 3.0.1 | bioconda | Model selection (ModelFinder) | ✓ INSTALLED |
| RAxML-NG | 1.2.2 | bioconda | Maximum likelihood tree search | ✓ INSTALLED |
| BEAST2 | 2.6.3 | bioconda | Bayesian inference | ✓ INSTALLED |
| BEAGLE | 4.0.1 | bioconda | BEAST2 likelihood acceleration | ✓ INSTALLED |
| MAFFT | 7.525 | bioconda | Multiple sequence alignment | ✓ INSTALLED |
| TrimAl | 1.5.1 | bioconda | Alignment trimming | ✓ INSTALLED |
| seqkit | — | bioconda | FASTA subsetting | ✓ INSTALLED |
| R | 4.5.2 | conda-forge | Statistical analysis | ✓ INSTALLED |
| r-caper | 1.0.4 | conda-forge | D-statistic (`phylo.d`) | ✓ INSTALLED |
| r-picante | 1.8.2 | conda-forge | `match.phylo.data()` | ✓ INSTALLED |
| r-phytools | 2.5-2 | conda-forge | Ancestral state reconstruction | ✓ INSTALLED |
| r-ape | 5.8-1 | conda-forge | Tree manipulation | ✓ INSTALLED |

> **Note on Java:** BEAST2 2.6.3 requires Java 8. Installing it downgraded `openjdk` from 25 to 8.0.144. IQ-TREE and RAxML-NG have been confirmed working after the downgrade.

---

## 11. Data Formats & File Specifications

Incorrect file formatting is the most common cause of silent failures — especially mismatched species names in the trait CSV.

### 11.1 Orthogroup Alignment FASTA

Files are named `OG*.trim.aln.fasta` and placed in `data/`.

```
>yHMPu5000035277_candida_pararugosa_201018.final
MSTEGKLVINGKPQVTIVGPTGLGKSTLL--EVNIVGVYEGDVSRQK...
>yHMPu5000035329_candida_saitoana_180604.haplomerger2.final
MSTEGKLVING-PQVTIVGPTGLGKSTLL--EVNIVGVYEGDVSRQK...
```

| Rule | Detail |
|------|--------|
| Header | Starts with `>` followed by species ID — no spaces, use underscores |
| Species names | Must match **exactly** the `Species_name_in_OGs` column in `mannose_traits.csv` — case sensitive |
| Sequence | Amino acid single-letter codes; gaps as `-`; all sequences same length (aligned) |
| Data type | Protein (AA) only — all scripts pass `-st AA` / `--data-type AA` |

### 11.2 Mannose Trait CSV

File location: `data/mannose_traits.csv`. Tab-separated, 3 columns, 267 data rows + 1 header.

```
Species	Species_name_in_OGs	Mannose_Class
Wickerhamiella_pararugosa	yHMPu5000035277_candida_pararugosa_201018.final	1
Candida_saitoana	yHMPu5000035329_candida_saitoana_180604.haplomerger2.final	1
Pichia_eremophila	yHMPu5000034641_pichia_eremophila_170307.haplomerger2.final	0
```

> **⚠ WARNING:** Species name mismatches are the most common failure point. `"Saccharomyces cerevisiae"` (space) vs `"Saccharomyces_cerevisiae"` (underscore) will silently drop that species. Always verify the `n_matched` column in `dstatistic_results.tsv` after running Phase 3.

### 11.3 jobs_with_models.txt (Phase 1 Output → Phase 2 Input)

Auto-generated by `collect_models.sh`. Tab-separated, no header. Do not edit manually.

```
OG0000007.trim.aln.fasta	Q.YEAST+F+R10
OG0000008.trim.aln.fasta	LG+F+R10
```

Verify before submitting Phase 2: `cat scripts/jobs_with_models.txt`

### 11.4 beast_jobs.txt (Phase 2b Input)

Generate after `generate_beast_xml.py` creates the XML files:

```bash
ls beast_xml/*.xml | xargs -n1 basename | sed 's/\.xml//' > scripts/beast_jobs.txt
```

One OG identifier per line, no extension, no path. Every name must have a corresponding `beast_xml/OG*.xml` file.

### 11.5 modelfinder_jobs.txt (Phase 1 Input)

Generate before submitting Phase 1:

```bash
ls data/OG*.trim.aln.fasta | xargs -n1 basename > scripts/modelfinder_jobs.txt
```

One alignment filename per line, basename only. To resubmit only failed OGs, manually edit the file to list only the failed filenames.

---

## 12. Functions Reference

### IQ-TREE ModelFinder

| Flag | Value | Explanation |
|------|-------|-------------|
| `-s` | `OG*.fasta` | Input alignment file |
| `-st AA` | `AA` | Declares amino acid data type (prevents misidentification as DNA) |
| `-m MF` | `MF` | Runs ModelFinder only — selects best model by BIC, does NOT build tree |
| `-T AUTO` | `AUTO` | Auto-detects optimal thread count |
| `--threads-max` | `4` | Caps threads at 4; matches `request_cpus` in the submit file |
| `--redo` | — | Overwrites existing checkpoints; prevents Condor retries from stalling |
| BIC | selection metric | Penalises model complexity; more conservative than AIC |

### RAxML-NG Tree Search

| Flag | Value | Explanation |
|------|-------|-------------|
| `--all` | — | Complete analysis: ML search + bootstrap + mapping in one step |
| `--data-type AA` | `AA` | Specifies amino acid data |
| `--model` | per-OG | Substitution model from `jobs_with_models.txt` (per-OG, not hardcoded) |
| `--tree pars{10},rand{10}` | 20 trees | 10 parsimony + 10 random starting trees; prevents local optima |
| `--bs-trees 100` | 100 | 100 bootstrap replicates; branch values ≥70% considered well-supported |
| `--seed 12345` | 12345 | Fixed seed for reproducibility |
| `--threads 10` | 10 | Must match `OMP_NUM_THREADS=10` in the Condor submit file |

### BEAST2 Bayesian Inference

| Parameter | Explanation |
|-----------|-------------|
| `UCRelaxedClockModel` | Relaxed lognormal clock — each branch gets its own rate from a lognormal distribution; more realistic for divergent yeasts |
| `BirthDeathGernhard08Model` | Birth-Death tree prior — models speciation and extinction; more realistic than Yule (no extinction) |
| `WAG site model` | Whelan & Goldman amino acid substitution matrix |
| `chainLength = 50,000,000` | 50M MCMC steps for sufficient posterior exploration |
| `logEvery = 5,000` | Records state every 5,000 steps → 10,000 samples in trace log |
| ESS > 200 | Effective Sample Size convergence threshold (ideally > 500) |
| `treeannotator -burnin 10` | Discards first 10% as burnin; summarizes remainder into MCC tree |

### R/ape, R/picante, R/caper

| Function | Explanation |
|----------|-------------|
| `read.tree()` | Reads Newick format tree into R `phylo` object (works for both RAxML and BEAST output) |
| `multi2di()` | Resolves polytomies to strictly binary tree (required by `phylo.d()`) |
| `chronos()` | Makes tree ultrametric (required by `phylo.d()` which assumes contemporary tips) |
| `match.phylo.data(tree, data)` | Links trait data to gene trees; returns only taxa present in **both**; unmatched taxa silently pruned |
| `comparative.data(phy, data, names.col)` | Builds comparative data object linking tree to trait data frame |
| `phylo.d(data, binvar, permut)` | Core D-statistic function; `permut=1000` random permutations for p-values |

---

## 13. D-statistic Interpretation

The D-statistic (Fritz & Purvis 2010) measures phylogenetic signal for binary traits, calculated using `caper::phylo.d()` in R — the same method used in the Y1000+ paper (Opulente et al. / Harrison et al.).

| D value | Biological interpretation |
|---------|--------------------------|
| D < 0 | Highly clustered — stronger signal than Brownian motion |
| D = 0 | Brownian motion — strong phylogenetic signal |
| 0 < D < 1 | Moderate phylogenetic signal |
| D = 1 | Random — no phylogenetic signal |
| D > 1 | Overdispersed — convergent evolution |

Two p-values are reported per analysis:

- **p(random)** — tests whether the trait is non-randomly distributed. Significant (< 0.05) means the trait has phylogenetic signal.
- **p(Brownian)** — tests whether the trait deviates from Brownian motion evolution. Significant means the trait is more or less clustered than expected under Brownian motion.

> **Note:** D-statistic on gene trees (this pipeline) answers a different question than D on the species tree. Gene tree D identifies which orthogroups are *associated* with mannose growth. Species tree D (planned future analysis) asks whether mannose growth is broadly conserved across all 267 yeast species.

---

## 14. Known Issues & Resolutions

| Issue | Resolution |
|-------|------------|
| Java downgrade | BEAST2 2.6.3 requires Java 8; installing it downgraded `openjdk` from 25 to 8.0.144. IQ-TREE and RAxML-NG confirmed working after downgrade. |
| Stray flag files | Early IQ-TREE test run created literal files named `-m`, `-pre`, `-st`, `-T` in project root due to a command formatting error. Removed with `rm -f -- '-m' '-pre' '-st' '-T'`. |
| Incomplete ModelFinder log | First ModelFinder run for OG0000007 was cut off mid-run (only 9/1,232 models tested). Incomplete logs deleted and run restarted with `nohup`. |
| Wrong `BASE_DIR` path | Scripts had hardcoded `/mnt/bigdata/linuxhome` path; actual path is `/home/glbrc.org/narhmadey`. Fixed with `sed -i` across all scripts. |
| `set -u` and `/etc/bashrc` | `BASHRCSOURCED: unbound variable` error from cluster `/etc/bashrc` when `set -u` is active. Fixed by wrapping conda activation with `set +u` / `set -u`. |
| `HISTTIMEFORMAT: unbound variable` | Interactive shell sets an unbound env var when using `set -u` in trim script. Fixed by adding `: "${HISTTIMEFORMAT:=}"` at the top of `trim_align.sh`. |
| OG0007179 removed | Outgroup (*Lipomyces starkeyi*) absent from this OG; fewer than 4 genes; cannot be properly rooted for D-statistic. Removed from Condor analysis. |

---

## 15. Current Status & Next Steps

| Phase | Status |
|-------|--------|
| Preprocessing (local) | ✅ COMPLETE — All OG FASTAs subset, aligned, and trimmed |
| Phase 1 — Model Selection | ✅ COMPLETE — All 8 OGs completed (OG0007179 removed) |
| Phase 2a — ML Trees | 🔄 IN PROGRESS — 8 RAxML-NG jobs running on cluster (cluster 105756) |
| Phase 2b — Bayesian Trees | ⏳ PENDING — MRCA calibration priors ready to inject into XML files |
| Phase 3 — D-statistic | ⏳ PENDING — `mannose_traits.csv` must be placed in `data/` before running |

### Quick Reference: Run in This Order

```bash
cd /home/glbrc.org/narhmadey/2.Y1000_267_phylogenetics/2_Plant_Path_563
conda activate phylo_env

# ── Phase 1: Submit ModelFinder jobs ──────────────────────────────
ls data/OG*.trim.aln.fasta | xargs -n1 basename > scripts/modelfinder_jobs.txt
condor_submit scripts/modelfinder_jobs.sub
condor_q                          # wait for all 9 jobs (~1.5 hrs)

# ── Phase 1: Collect models after jobs finish ──────────────────────
bash scripts/collect_models.sh
cat scripts/jobs_with_models.txt  # verify all models present

# ── Phase 2a + 2b: Submit simultaneously ──────────────────────────
condor_submit scripts/raxml_jobs.sub
python3 scripts/generate_beast_xml.py
ls beast_xml/*.xml | xargs -n1 basename | sed 's/\.xml//' > scripts/beast_jobs.txt
condor_submit scripts/beast_jobs.sub
condor_q                          # monitor both sets of jobs

# ── Phase 2b: Summarize Bayesian trees after BEAST2 jobs finish ───
for f in results/beast/*.trees; do
    og=$(basename $f .trees)
    treeannotator -burnin 10 $f results/beast/${og}.beast.tree
done

# ── Phase 3: D-statistic ──────────────────────────────────────────
# Place mannose_traits.csv in data/ first, then:
Rscript scripts/run_dstatistic.R \
    --trait data/mannose_traits.csv \
    --treedir results/ \
    --output results/dstatistic_results.tsv
```

### Planned Next Analyses

1. Per-site conservation scoring (Shannon entropy, Jensen-Shannon divergence)
2. Mannose grower vs. non-grower enrichment tests per OG
3. Gene family presence/absence modeling
4. Mannose metabolism pathway mapping
5. ASTRAL species-tree inference from ML gene trees
6. Ancestral state reconstruction for mannose growth trait
7. Genome alignment with pygenomeviz/progressiveMauve (dataset recovery in progress)

---

*GLBRC Cluster: scarcity-ap-1 | User: narhmadey | Conda env: phylo_env | Last updated: March 2026*
