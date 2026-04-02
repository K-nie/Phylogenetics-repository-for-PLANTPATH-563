# Phylogenetics-repository-for-PLANTPATH-563 (This repo was created for the purposes for this class). 
Mannose-associated conserved genes in 267 yeast species. Phylogenomic pipeline to identify and analyze gene families associated with mannose utilization across 267 yeast species using orthogroups, 
multiple sequence alignment, trimming, and gene tree inference.

## Project overview
This project analyzes 267 focal yeast species drawn from the Y1000+ dataset to identify genes conserved in mannose-utilizing species compared with non-utilizers.
The workflow:
1. Quality control of input data
2. Build comparable gene sets via orthogroups
3. Subset orthogroups to 267 species
4. Perform multiple sequence alignment (MSA) per orthogroup
5. Trim low-quality alignment regions
6. Infer gene trees
7. Quantify conservation and test phenotype associations

## Core biological logic
Machine learning models identified orthogroups (OGs) associated with mannose phenotype. Each orthogroup represents a gene family. Each species may contribute 0, 1, or multiple genes per OG.
Therefore, MSA is performed per orthogroup, not per species.

## Repository structure
data/
├── raw/
│   └── y1000p_orthofinder/
│       └── Orthogroup_Sequences/        # Full OG FASTAs (1154 species)
│
└── processed/
    └── 03_ml_orthogroups/
        ├── 01_og_lists/
        │   └── orthogroup_list.txt      # OGs from ML models
        │
        ├── 05_og_fastas_267/            # OG FASTAs subset to 267 species
        ├── 06_alignments/               # MAFFT outputs
        ├── 07_trimmed_alignments/       # TrimAl outputs
        ├── 08_gene_trees/               # IQ-TREE outputs
        └── logs/                        # pipeline logs

## Environment setup
Create the reproducible phylogenomics environment:

conda create -n phylo_env -c bioconda -c conda-forge \
    iqtree mafft trimal fasttree seqkit
conda activate phylo_env

## Input data
### Required inputs
File	                        Description
Orthogroup_Sequences/*.fasta	Orthofinder OG sequences (1154 species)
species_267.txt	                List of focal species
orthogroup_list.txt	            OGs selected by ML models

# Step 1 — Subset orthogroups to 267 species
Each OG FASTA is filtered to retain only sequences from the focal species.
Script: og_fastas_267.sh

# Step 2 — Multiple sequence alignment (MAFFT)
Each orthogroup is aligned independently.
## Strategy used
OG size	        MAFFT mode	    Reason
≤200 sequences	L-INS-i	        high accuracy
>200 sequences	--auto	        scalable

# Step 3 — Alignment trimming (TrimAl)
Poorly aligned regions are removed prior to phylogenetic inference.
Script: trim_align.sh

# Step 4 — Gene tree inference (IQ-TREE)
After trimming, gene trees are inferred.
Example command
iqtree3 -s alignment.trim.aln.fasta \
        -m MFP \
        -B 1000 \
        -T AUTO \
        -pre output_prefix

## Parameter rationale
Parameter	    Purpose
-m MFP	        automatic model selection
-B 1000	        ultrafast bootstrap support
-T AUTO	        automatic threading

## Next analyses (planned)
1. Per-site conservation scoring
2. Grower vs non-grower enrichment tests
3. Gene family presence/absence modeling
4. Mannose metabolism pathway mapping
5. Possible species-tree inference
