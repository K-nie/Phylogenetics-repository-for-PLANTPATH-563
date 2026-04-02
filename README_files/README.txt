  
Saccharomycotina Yeast Phylogenetics Pipeline

  README — Analysis Documentation
  Project: 2_Plant_Path_563  |  Author: Benjamin Narh-Madey  |  Environment: phylo_env  |  Date: March 2026
1. Project Overview
This pipeline performs phylogenetic analysis on aligned orthogroup protein sequences from the Y1000+ Saccharomycotina yeast dataset. The analysis runs on the GLBRC HT Condor cluster and consists of three sequential phases.

The dataset currently contains 9 orthogroup alignments (267 species) located in the data/ directory. A test run was completed on OG0000007 to validate the pipeline before full batch submission.

1.1 Research Questions
	•	Which substitution model best fits each orthogroup alignment?
	•	What is the maximum likelihood gene tree topology for each orthogroup?
	•	What is the Bayesian posterior gene tree for each orthogroup?
	•	Does mannose growth ability show phylogenetic signal across yeast gene trees?

1.2 Dataset
	•	9 orthogroup protein alignments (OG*.trim.aln.fasta) in data/
	•	267 yeast species from the Y1000+ Saccharomycotina dataset
	•	Binary mannose trait data (1 = growth, 0 = no growth) for D-statistic analysis
	•	Source data in 00_downloads/y1000p_orthofinder.tar.gz and y1000p_pep_files.tar.gz

2. Pipeline Overview
The pipeline is divided into three phases. Phases 2a and 2b run in parallel on Condor after Phase 1 completes.

Phase 1 — Model Selection  (submit node, serial, ~13 hrs for 9 OGs)
run_modelfinder.sh  →  jobs_with_models.txt

Phase 2a — ML Trees  (Condor, 9 parallel jobs, ~2 hrs each)
raxml_jobs.sub  →  run_raxml_job.sh  →  results/OG*.raxml.support

Phase 2b — Bayesian Trees  (Condor, 9 parallel jobs, ~24 hrs each)
beast_jobs.sub  →  run_beast_job.sh  →  results/beast/OG*.beast.tree

Phase 3 — D-statistic  (submit node, after Phase 2 complete)
run_dstatistic.R  →  results/dstatistic_results.tsv

3. Step-by-Step Instructions
Phase 1 — Model Selection
Run on the submit node (scarcity-ap-1). Use nohup so it survives SSH disconnections.

# Activate environment
conda activate phylo_env

# Launch model selection in background
nohup bash scripts/run_modelfinder.sh > logs/modelfinder.out 2>&1 &
echo "PID: $!"

# Monitor progress
tail -f logs/modelfinder.out

# Check which models have been written so far
cat scripts/jobs_with_models.txt

When complete, jobs_with_models.txt will contain one line per orthogroup:
OG0000007.trim.aln.fasta    Q.YEAST+F+R10
OG0000008.trim.aln.fasta    LG+R8
...

Phase 2a — Maximum Likelihood Trees (RAxML-NG)
Submit all 9 RAxML jobs simultaneously after Phase 1 completes.

# Review job list first
cat scripts/jobs_with_models.txt

# Submit all 9 ML jobs at once
condor_submit scripts/raxml_jobs.sub

# Monitor jobs
condor_q

Each job runs on 10 CPUs with 8GB RAM. Settings: 10 parsimony + 10 random starting trees, 100 bootstrap replicates. Results written to results/.

Phase 2b — Bayesian Inference (BEAST2)
Runs in parallel with Phase 2a. Settings: relaxed lognormal clock, Birth-Death tree prior, 50 million MCMC steps.

# Step 1: Generate one XML config per orthogroup
python3 scripts/generate_beast_xml.py

# Step 2: Generate job list
ls beast_xml/*.xml | xargs -n1 basename | sed 's/\.xml//' > scripts/beast_jobs.txt

# Step 3: Submit all 9 BEAST2 jobs
condor_submit scripts/beast_jobs.sub

# After jobs complete — summarize posterior trees (10% burnin)
for f in results/beast/*.trees; do
  og=$(basename $f .trees)
  treeannotator -burnin 10 $f results/beast/${og}.beast.tree
done

IMPORTANT: Before using Bayesian trees, check convergence in Tracer. All ESS values should be > 200. If not, the chain length needs to be increased.

Phase 3 — D-statistic Analysis
Run after Phase 2 trees are available. Calculates phylogenetic signal for mannose growth trait across gene trees.

Rscript scripts/run_dstatistic.R \
  --trait data/mannose_traits.csv \
  --treedir results/ \
  --output results/dstatistic_results.tsv

The mannose_traits.csv file must have two columns: species (matching tree tip labels exactly) and mannose (0 or 1).

4. Scripts Reference
Script
Phase
Description
run_modelfinder.sh
Phase 1
Runs IQ-TREE ModelFinder serially on all orthogroups. Generates jobs_with_models.txt
raxml_jobs.sub
Phase 2a
Condor submit file for RAxML-NG ML tree search jobs
run_raxml_job.sh
Phase 2a
Per-job script: unpacks env, runs RAxML-NG, copies results
beast_jobs.sub
Phase 2b
Condor submit file for BEAST2 Bayesian inference jobs
run_beast_job.sh
Phase 2b
Per-job script: unpacks env, runs BEAST2, copies results
generate_beast_xml.py
Phase 2b
Generates one BEAST2 XML config per orthogroup alignment
run_dstatistic.R
Phase 3
Calculates D-statistic for mannose trait on gene trees

5. Output Files
Directory
File Pattern
Phase
Description
iqtree_logs/
OG*.trim.aln_model.log
Phase 1
IQ-TREE ModelFinder full log per orthogroup
scripts/
jobs_with_models.txt
Phase 1
Two-column file: filename + best model (input to Phase 2)
results/
OG*.raxml.bestTree
Phase 2a
Best ML tree topology (no support values)
results/
OG*.raxml.support
Phase 2a
Best ML tree with bootstrap support values — MAIN RESULT
results/
OG*.raxml.bestModel
Phase 2a
Best-fit substitution model used by RAxML
results/
OG*.raxml.bootstraps
Phase 2a
All 100 individual bootstrap trees
results/beast/
OG*.trace.log
Phase 2b
BEAST2 MCMC trace — check in Tracer (ESS > 200)
results/beast/
OG*.trees
Phase 2b
Posterior tree distribution from BEAST2
results/beast/
OG*.beast.tree
Phase 2b
Summarized Bayesian tree (after TreeAnnotator)
results/
dstatistic_results.tsv
Phase 3
D-statistic per orthogroup: D value, p-values, interpretation
logs/
*.out / *.err
All
Condor stdout/stderr logs per job

6. Tools & Environment
All tools are installed in the phylo_env conda environment at:
/home/glbrc.org/narhmadey/.conda/envs/phylo_env

The environment is also packaged as phylo_env.tar.gz for deployment to Condor execute nodes.

Tool
Version
Channel
Purpose
Status
IQ-TREE
3.0.1
bioconda
Model selection (ModelFinder)
✓ INSTALLED
RAxML-NG
1.2.2
bioconda
Maximum likelihood tree search
✓ INSTALLED
BEAST2
2.6.3
bioconda
Bayesian inference
✓ INSTALLED
BEAGLE
4.0.1
bioconda
BEAST2 likelihood acceleration
✓ INSTALLED
R
4.5.2
conda-forge
Statistical analysis
✓ INSTALLED
r-caper
1.0.4
conda-forge
D-statistic (phylo.d)
✓ INSTALLED
r-picante
1.8.2
conda-forge
match.phylo.data()
✓ INSTALLED
r-phytools
2.5-2
conda-forge
Ancestral state reconstruction
✓ INSTALLED
r-ape
5.8-1
conda-forge
Tree manipulation
✓ INSTALLED
MAFFT
7.525
bioconda
Multiple sequence alignment
✓ INSTALLED
TrimAl
1.5.1
bioconda
Alignment trimming
✓ INSTALLED

6.1 Activate Environment
conda activate phylo_env

6.2 Repackage Environment (after installing new tools)
conda pack -n phylo_env -o phylo_env.tar.gz

7. D-statistic Interpretation
The D-statistic (Fritz & Purvis 2010) measures phylogenetic signal for binary traits. It is calculated using caper::phylo.d() in R, the same method used in the Y1000+ paper (Opulente et al. / Harrison et al.).

D value
Biological interpretation
D < 0
Highly clustered — stronger than Brownian motion
D = 0
Brownian motion — strong phylogenetic signal
0 < D < 1
Moderate phylogenetic signal
D = 1
Random — no phylogenetic signal
D > 1
Overdispersed — convergent evolution

Two p-values are reported per analysis:
	•	p(random) — tests whether the trait is non-randomly distributed. Significant (< 0.05) means the trait has phylogenetic signal.
	•	p(Brownian) — tests whether the trait deviates from Brownian motion evolution. Significant means the trait is more or less clustered than expected under Brownian motion.

Note: D-statistic on gene trees (this pipeline) answers a different question than D on the species tree. Gene tree D identifies which orthogroups are associated with mannose growth. Species tree D (future analysis) asks whether mannose growth is broadly conserved across all 267 yeast species.

8. Directory Structure
2_Plant_Path_563/
├── data/                    # Input alignments (OG*.trim.aln.fasta)
├── beast_xml/               # BEAST2 XML configs (one per OG)
├── iqtree_logs/             # IQ-TREE ModelFinder output
├── results/                 # RAxML-NG output trees
│   └── beast/               # BEAST2 output trees
├── logs/                    # Condor + runtime logs
├── scripts/                 # All pipeline scripts
│   ├── run_modelfinder.sh
│   ├── raxml_jobs.sub
│   ├── run_raxml_job.sh
│   ├── beast_jobs.sub
│   ├── run_beast_job.sh
│   ├── generate_beast_xml.py
│   ├── run_dstatistic.R
│   ├── jobs_with_models.txt # Generated by Phase 1
│   └── beast_jobs.txt       # Generated before Phase 2b
├── bin/                     # raxml-ng binary
├── 00_downloads/            # Raw downloaded data
└── phylo_env.tar.gz         # Packaged conda environment

9. Data Formats & File Specifications
This section describes the exact format required for each input and output file in the pipeline. Incorrect formatting is the most common cause of silent failures — especially mismatched species names in the trait CSV.

9.1 Orthogroup Alignment FASTA (Input to Phase 1 & 2)
Each orthogroup alignment file must be a trimmed, multiple sequence alignment in FASTA format. Files are named OG*.trim.aln.fasta and placed in the data/ directory.

Required filename pattern:
data/OG0000007.trim.aln.fasta
data/OG0000008.trim.aln.fasta
data/OG*.trim.aln.fasta   ← wildcard used by all scripts

File format — each sequence on two lines:
>Saccharomyces_cerevisiae
MSTEGKLVINGKPQVTIVGPTGLGKSTLLNRLLGRDSPTEG--EVNIVGVYEGDVSRQK
>Candida_albicans
MSTEGKLVING-PQVTIVGPTGLGKSTLLNRLLGRDSPTEG--EVNIVGVYEGDVSRQK
>Kluyveromyces_marxianus
MSTEGKLVINGKPQVTIVGPTGLGKSTLLNRLLGRDSPTE---EVNIVGVYEGDVSRQK

Rule
Detail
Header line
Starts with > followed immediately by species name — no spaces in name, use underscores
Species names
Must match EXACTLY the names in mannose_traits.csv — case sensitive, underscores not spaces
Sequence
Amino acid single-letter codes. Gaps represented as - (dash). All sequences must be same length (aligned)
Data type
Protein (AA) only — scripts pass -st AA / --data-type AA to all tools
Trimming
Alignments should already be trimmed with TrimAl (installed). The .trim. in the filename indicates this.

9.2 Mannose Trait CSV (Input to Phase 3 D-statistic)
The mannose trait file is a comma-separated CSV with two columns. It must be saved as data/mannose_traits.csv before running Phase 3.

Required file location:
data/mannose_traits.csv

File format:
species,mannose
Saccharomyces_cerevisiae,1
Candida_albicans,0
Kluyveromyces_marxianus,1
Yarrowia_lipolytica,0
Pichia_kudriavzevii,1
...267 rows total

Rule
Detail
Column 1: species
Species name — must match EXACTLY the tip labels in gene trees produced by RAxML/BEAST. Use underscores, not spaces.
Column 2: mannose
Binary trait: 1 = growth on mannose, 0 = no growth. No other values allowed.
Header row
Must be present: species,mannose — the R script auto-detects these column names
Missing data
Do not include rows with missing mannose values. Remove those species entirely from the CSV.
Row count
267 rows of data (plus 1 header row = 268 lines total). Check with: wc -l data/mannose_traits.csv
Name matching
The R script uses picante::match.phylo.data() to match names. Species in CSV not found in tree are silently dropped. Always verify match count in output.

WARNING: Species name mismatches are the most common failure point. A name like "Saccharomyces cerevisiae" (space) vs "Saccharomyces_cerevisiae" (underscore) will silently drop that species. Always check the n_matched column in dstatistic_results.tsv after running Phase 3.

9.3 jobs_with_models.txt (Phase 1 Output → Phase 2 Input)
This file is automatically generated by run_modelfinder.sh at the end of Phase 1. It feeds directly into raxml_jobs.sub as the job queue. You should not edit it manually — but you can inspect it to verify model selection results before submitting Phase 2.

File location:
scripts/jobs_with_models.txt

File format — tab-separated, no header:
# Column 1: alignment filename    Column 2: best-fit model (tab-separated)
OG0000007.trim.aln.fasta	Q.YEAST+F+R10
OG0000008.trim.aln.fasta	LG+R8
OG0000009.trim.aln.fasta	WAG+F+R6
OG0000020.trim.aln.fasta	Q.PFAM+G4
OG0000029.trim.aln.fasta	Q.YEAST+R10
OG0000032.trim.aln.fasta	LG+F+R7
OG0000044.trim.aln.fasta	WAG+G4
OG0000067.trim.aln.fasta	Q.INSECT+F+R8
OG0007179.trim.aln.fasta	LG+R6

Rule
Detail
Separator
Tab-separated (\t) — NOT comma. raxml_jobs.sub reads it with: queue og_file, best_model from jobs_with_models.txt
No header row
File contains only data lines — no column names. Condor reads it directly.
Line count
Must equal number of alignment files in data/. Check with: wc -l scripts/jobs_with_models.txt
Model format
Exact IQ-TREE model string passed directly to RAxML-NG --model flag. Do not modify.
Verify before submit
Always cat scripts/jobs_with_models.txt before running condor_submit to confirm all 9 OGs have models.

9.4 beast_jobs.txt (Phase 2b Input)
This file lists the orthogroup names to submit as BEAST2 jobs. Unlike jobs_with_models.txt, it is generated manually with a one-liner after generate_beast_xml.py creates the XML files.

File location:
scripts/beast_jobs.txt

Generate it with this command:
ls beast_xml/*.xml | xargs -n1 basename | sed 's/\.xml//' > scripts/beast_jobs.txt

File format — one orthogroup name per line, no extension, no path:
# One orthogroup name per line — no .xml extension, no path
OG0000007
OG0000008
OG0000009
OG0000020
OG0000029
OG0000032
OG0000044
OG0000067
OG0007179

Rule
Detail
One name per line
Just the OG identifier — no file extension, no path, no model. beast_jobs.sub passes this as $1 to run_beast_job.sh.
No header row
File contains only OG names. Condor reads each line as the og_name variable.
Matches XML files
Every name in beast_jobs.txt must have a corresponding XML file in beast_xml/. E.g. OG0000007 requires beast_xml/OG0000007.xml
Partial runs
To resubmit only failed jobs, manually edit beast_jobs.txt to list only the failed OG names before rerunning condor_submit.

10. Known Issues & Resolutions
Issue
Resolution
Java downgrade
BEAST2 2.6.3 requires Java 8; installing it downgraded openjdk from 25 to 8.0.144. IQ-TREE and RAxML-NG confirmed working after downgrade.
Stray flag files
Early IQ-TREE test run created literal files named -m, -pre, -st, -T in project root due to command formatting issue. Removed with: rm -f -- '-m' '-pre' '-st' '-T'
Incomplete ModelFinder log
First ModelFinder run for OG0000007 was cut off mid-run (only 9/1232 models tested). Deleted incomplete logs and restarted with nohup.
Wrong BASE_DIR path
Scripts had hardcoded /mnt/bigdata/linuxhome path; actual path is /home/glbrc.org/narhmadey. Fixed with sed -i across all scripts.
set -u and /etc/bashrc
BASHRCSOURCED unbound variable error from cluster /etc/bashrc when set -u active. Fixed by wrapping conda activation with set +u / set -u.

11. Current Status
Phase 1 — Model Selection: IN PROGRESS — ModelFinder running via nohup on submit node
Phase 2a — ML Trees: PENDING — waiting for Phase 1 to complete
Phase 2b — Bayesian Trees: PENDING — generate_beast_xml.py ready to run
Phase 3 — D-statistic: PENDING — mannose_traits.csv needs to be placed in data/

Test run completed successfully on OG0000007 using IQ-TREE alone. Best-fit model: Q.YEAST+F+R10 (BIC). Full pipeline scripts written and path issues resolved.

GLBRC Cluster: scarcity-ap-1  |  User: narhmadey  |  Conda env: phylo_env  |  Last updated: March 2026
