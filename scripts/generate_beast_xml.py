#!/usr/bin/env python3
"""
generate_beast_xml.py
Generates BEAST2 XML configuration files per orthogroup alignment.
UPGRADED: Dynamically matches IQ-TREE models and generates 2 runs per orthogroup.

USAGE:
    python3 scripts/generate_beast_xml.py

OUTPUT:
    beast_xml/OG0000007_run1.xml
    beast_xml/OG0000007_run2.xml
    ...two XMLs per orthogroup alignment in data/

REQUIRES:
    - Biopython
    - Alignment files in data/OG*.trim.aln.fasta
    - jobs_with_models.txt file with IQ-TREE best models
Author: Benjamin Narh-Madey
"""

import os
from pathlib import Path
from Bio import AlignIO

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR    = Path("/home/glbrc.org/narhmadey/2.Y1000_267_phylogenetics/2_Plant_Path_563")
ALIGN_DIR   = BASE_DIR / "data"
XML_DIR     = BASE_DIR / "beast_xml"
MODEL_FILE  = BASE_DIR / "scripts" / "jobs_with_models.txt"  # per-OG IQ-TREE models
IQTREE_LOGS = BASE_DIR / "iqtree_logs"                       # IQ-TREE Phase 1 output
XML_DIR.mkdir(parents=True, exist_ok=True)

# ── BEAST2 settings ───────────────────────────────────────────────────────────
CHAIN_LENGTH  = 50_000_000
LOG_EVERY     = 5_000        # log every 5000 steps → 10000 samples total
TRACE_EVERY   = 5_000
TREE_EVERY    = 5_000


def read_alignment(fasta_path):
    """Read alignment and return list of (taxon, sequence) tuples."""
    aln = AlignIO.read(fasta_path, "fasta")
    return [(rec.id, str(rec.seq)) for rec in aln]


def read_models(model_filepath):
    """Parse jobs_with_models.txt to map orthogroups to their best models."""
    og_models = {}
    if not model_filepath.exists():
        print(f"WARNING: Model file {model_filepath} not found. Defaulting all to WAG.")
        return og_models
        
    with open(model_filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    # Clean filename to get just the OG name
                    og = parts[0].replace(".trim.aln.fasta", "").replace(".fasta", "")
                    
                    # IQ-TREE models look like "LG+F+G4" or "Q.YEAST+I+R4"
                    # We just need the base matrix name (before the first '+')
                    full_model = parts[1]
                    base_model = full_model.split('+')[0]
                    og_models[og] = base_model
    return og_models


def make_sequence_block(taxa_seqs):
    """Generate BEAST2 sequence XML block."""
    lines = []
    for taxon, seq in taxa_seqs:
        lines.append(f'        <sequence id="seq_{taxon}" taxon="{taxon}" '
                     f'totalcount="20" value="{seq}"/>')
    return "\n".join(lines)


def make_taxon_block(taxa_seqs):
    """Generate taxon set XML."""
    lines = []
    for taxon, _ in taxa_seqs:
        lines.append(f'                    <taxon id="{taxon}" spec="Taxon"/>')
    return "\n".join(lines)


def read_iqtree_starting_tree(og_name):
    """
    Read the IQ-TREE ML tree produced in Phase 1 (-m MFP) for use as
    a BEAST2 starting tree. This improves MCMC convergence by starting
    from a reasonable topology rather than a random tree.
    Returns the Newick string if found, or None to fall back to RandomTree.
    """
    treefile = IQTREE_LOGS / f"{og_name}.trim.aln_model.treefile"
    if treefile.exists():
        with open(treefile) as f:
            newick = f.read().strip()
        if newick:
            print(f"    Using IQ-TREE starting tree from: {treefile.name}")
            return newick
    print(f"    No IQ-TREE treefile found for {og_name} — using RandomTree initialiser")
    return None


def make_init_block(og_name, taxa_seqs, starting_newick):
    """
    Generate the BEAST2 tree initialisation block.
    - If an IQ-TREE ML tree is available, use TreeParser with the Newick string.
      This seeds BEAST2 from a good starting topology, reducing burn-in.
    - Otherwise fall back to a random coalescent tree (ConstantPopulation model).
    """
    if starting_newick:
        return f"""        <init estimate="false" id="RandomTree.t:{og_name}"
              initial="@Tree.t:{og_name}" spec="beast.util.TreeParser"
              taxa="@{og_name}" IsLabelledNewick="true"
              newick="{starting_newick}"/>"""
    else:
        return f"""        <init estimate="false" id="RandomTree.t:{og_name}"
              initial="@Tree.t:{og_name}" spec="beast.evolution.tree.RandomTree"
              taxa="@{og_name}">
            <populationModel spec="beast.evolution.tree.coalescent.ConstantPopulation"
                             popSize="0.1"/>
        </init>"""


def generate_xml(og_name, taxa_seqs, base_model, run_id, out_path, starting_newick=None):
    """Generate complete BEAST2 XML for one orthogroup."""
    n_taxa      = len(taxa_seqs)
    seq_block   = make_sequence_block(taxa_seqs)
    taxon_block = make_taxon_block(taxa_seqs)
    init_block  = make_init_block(og_name, taxa_seqs, starting_newick)

    xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast beautitemplate="Standard" beautistatus="" namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" required="" version="2.6">

    <data id="{og_name}" dataType="aminoacid" name="alignment">
{seq_block}
    </data>

    <map name="Uniform">beast.math.distributions.Uniform</map>
    <map name="Exponential">beast.math.distributions.Exponential</map>
    <map name="LogNormal">beast.math.distributions.LogNormalDistributionModel</map>
    <map name="Normal">beast.math.distributions.Normal</map>
    <map name="Beta">beast.math.distributions.Beta</map>
    <map name="Gamma">beast.math.distributions.Gamma</map>
    <map name="LaplaceDistribution">beast.math.distributions.LaplaceDistribution</map>
    <map name="prior">beast.math.distributions.Prior</map>
    <map name="InverseGamma">beast.math.distributions.InverseGamma</map>
    <map name="OneOnX">beast.math.distributions.OneOnX</map>

    <run chainLength="{CHAIN_LENGTH}" id="mcmc" spec="MCMC">

        <state id="state" storeEvery="5000">

            <tree id="Tree.t:{og_name}" name="stateNode">
                <taxonset id="TaxonSet.{og_name}" spec="TaxonSet">
{taxon_block}
                </taxonset>
            </tree>

            <parameter id="BDBirthRate.t:{og_name}"
                       lower="0.0" name="stateNode" upper="10000.0">1.0</parameter>
            <parameter id="BDDeathRate.t:{og_name}"
                       lower="0.0" name="stateNode" upper="1.0">0.5</parameter>

            <parameter id="ucldMean.c:{og_name}"
                       lower="0.0" name="stateNode">0.001</parameter>
            <parameter id="ucldSdev.c:{og_name}"
                       lower="0.0" name="stateNode" upper="10.0">0.1</parameter>
            <parameter id="rateCategories.c:{og_name}"
                       dimension="{2 * n_taxa - 2}" name="stateNode">1</parameter>

            <parameter id="freqParameter.s:{og_name}"
                       dimension="20" lower="0.0" name="stateNode" upper="1.0">0.05</parameter>

            <parameter id="gammaShape.s:{og_name}"
                       lower="0.0" name="stateNode">0.5</parameter>
        </state>

{init_block}

        <distribution id="posterior" spec="util.CompoundDistribution">

            <distribution id="prior" spec="util.CompoundDistribution">

                <distribution
                    id="BirthDeath.t:{og_name}"
                    spec="beast.evolution.speciation.BirthDeathGernhard08Model"
                    birthDiffRate="@BDBirthRate.t:{og_name}"
                    relativeDeathRate="@BDDeathRate.t:{og_name}"
                    tree="@Tree.t:{og_name}"/>

                <prior id="BDBirthRatePrior.t:{og_name}"
                       name="distribution" x="@BDBirthRate.t:{og_name}">
                    <LogNormal M="1.0" S="1.25" meanInRealSpace="false"
                               name="distr"/>
                </prior>

                <prior id="BDDeathRatePrior.t:{og_name}"
                       name="distribution" x="@BDDeathRate.t:{og_name}">
                    <Beta alpha="1.0" beta="1.0" name="distr"/>
                </prior>

                <prior id="ucldMeanPrior.c:{og_name}"
                       name="distribution" x="@ucldMean.c:{og_name}">
                    <LogNormal M="-7.0" S="1.0" meanInRealSpace="false"
                               name="distr"/>
                </prior>

                <prior id="ucldSdevPrior.c:{og_name}"
                       name="distribution" x="@ucldSdev.c:{og_name}">
                    <Exponential mean="0.3333" name="distr"/>
                </prior>

                <prior id="GammaShapePrior.s:{og_name}"
                       name="distribution" x="@gammaShape.s:{og_name}">
                    <Exponential mean="0.1" name="distr"/>
                </prior>

            </distribution>

            <distribution id="likelihood" spec="util.CompoundDistribution"
                          useThreads="true">
                <distribution id="treeLikelihood.{og_name}"
                              spec="ThreadedTreeLikelihood"
                              data="@{og_name}" tree="@Tree.t:{og_name}">

                    <siteModel id="SiteModel.s:{og_name}"
                               spec="SiteModel"
                               gammaCategoryCount="4"
                               shape="@gammaShape.s:{og_name}">
                        <substModel id="{base_model}.s:{og_name}"
                                    spec="{base_model}"
                                    frequencies="@freqParameter.s:{og_name}"/>
                    </siteModel>

                    <branchRateModel id="RelaxedClock.c:{og_name}"
                                     spec="beast.evolution.branchratemodel.UCRelaxedClockModel"
                                     rateCategories="@rateCategories.c:{og_name}"
                                     tree="@Tree.t:{og_name}">
                        <LogNormal S="@ucldSdev.c:{og_name}"
                                   M="@ucldMean.c:{og_name}"
                                   meanInRealSpace="true" name="distr"/>
                    </branchRateModel>

                </distribution>
            </distribution>
        </distribution>

        <operator id="CoalescentExchangeNarrow.t:{og_name}"
                  spec="Exchange" tree="@Tree.t:{og_name}" weight="15.0"/>
        <operator id="CoalescentExchangeWide.t:{og_name}"
                  spec="Exchange" isNarrow="false"
                  tree="@Tree.t:{og_name}" weight="3.0"/>
        <operator id="SubtreeSlide.t:{og_name}"
                  spec="SubtreeSlide" tree="@Tree.t:{og_name}" weight="15.0"/>
        <operator id="TreeUniform.t:{og_name}"
                  spec="Uniform" tree="@Tree.t:{og_name}" weight="30.0"/>
        <operator id="TreeRootScaler.t:{og_name}"
                  spec="ScaleOperator" rootOnly="true"
                  scaleFactor="0.95" tree="@Tree.t:{og_name}" weight="3.0"/>
        <operator id="TreeScaler.t:{og_name}"
                  spec="ScaleOperator" scaleFactor="0.95"
                  tree="@Tree.t:{og_name}" weight="3.0"/>
        <operator id="WilsonBalding.t:{og_name}"
                  spec="WilsonBalding" tree="@Tree.t:{og_name}" weight="3.0"/>

        <operator id="BDBirthRateScaler.t:{og_name}"
                  spec="ScaleOperator" parameter="@BDBirthRate.t:{og_name}"
                  scaleFactor="0.75" weight="3.0"/>
        <operator id="BDDeathRateScaler.t:{og_name}"
                  spec="ScaleOperator" parameter="@BDDeathRate.t:{og_name}"
                  scaleFactor="0.75" weight="3.0"/>

        <operator id="ucldMeanScaler.c:{og_name}"
                  spec="ScaleOperator" parameter="@ucldMean.c:{og_name}"
                  scaleFactor="0.75" weight="3.0"/>
        <operator id="ucldSdevScaler.c:{og_name}"
                  spec="ScaleOperator" parameter="@ucldSdev.c:{og_name}"
                  scaleFactor="0.75" weight="3.0"/>
        <operator id="rateCategoriesRandomWalk.c:{og_name}"
                  spec="IntRandomWalkOperator"
                  parameter="@rateCategories.c:{og_name}"
                  windowSize="1" weight="10.0"/>
        <operator id="rateCategoriesSwapOperator.c:{og_name}"
                  spec="SwapOperator"
                  intparameter="@rateCategories.c:{og_name}" weight="10.0"/>
        <operator id="rateCategoriesUniform.c:{og_name}"
                  spec="UniformOperator"
                  parameter="@rateCategories.c:{og_name}" weight="10.0"/>

        <operator id="GammaShapeScaler.s:{og_name}"
                  spec="ScaleOperator" parameter="@gammaShape.s:{og_name}"
                  scaleFactor="0.75" weight="0.1"/>

        <logger id="screenlog" logEvery="{LOG_EVERY}" spec="Logger">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log id="ESS.0" spec="util.ESS" arg="@posterior" name="log"/>
        </logger>

        <logger fileName="{og_name}_{run_id}.trace.log" id="tracelog"
                logEvery="{TRACE_EVERY}" model="@posterior"
                sanitiseHeaders="true" sort="smart" spec="Logger">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood.{og_name}"/>
            <log idref="BirthDeath.t:{og_name}"/>
            <parameter idref="BDBirthRate.t:{og_name}" name="log"/>
            <parameter idref="BDDeathRate.t:{og_name}" name="log"/>
            <parameter idref="ucldMean.c:{og_name}" name="log"/>
            <parameter idref="ucldSdev.c:{og_name}" name="log"/>
            <parameter idref="gammaShape.s:{og_name}" name="log"/>
            <log id="rate.c:{og_name}"
                 spec="beast.evolution.branchratemodel.RateStatistic"
                 branchRateModel="@RelaxedClock.c:{og_name}"
                 tree="@Tree.t:{og_name}"/>
        </logger>

        <logger fileName="{og_name}_{run_id}.trees" id="treelog"
                logEvery="{TREE_EVERY}" mode="tree" spec="Logger">
            <log id="TreeWithMetaDataLogger.t:{og_name}"
                 spec="beast.evolution.tree.TreeWithMetaDataLogger"
                 branchRateModel="@RelaxedClock.c:{og_name}"
                 tree="@Tree.t:{og_name}"/>
        </logger>

    </run>
</beast>
"""
    with open(out_path, "w") as fh:
        fh.write(xml)


def main():
    alignments = sorted(ALIGN_DIR.glob("OG*.trim.aln.fasta"))
    if not alignments:
        print(f"ERROR: No alignments found in {ALIGN_DIR}")
        raise SystemExit(1)

    # 1. Load the specific models found by IQ-TREE
    og_models = read_models(MODEL_FILE)

    print(f"Generating BEAST2 XML files ({len(alignments)} orthogroups x 2 runs each)")
    print(f"  Settings: Relaxed lognormal clock | Birth-Death prior | {CHAIN_LENGTH:,} MCMC steps")
    print(f"  Output:   {XML_DIR}/")
    print()

    for aln_path in alignments:
        og_name = aln_path.name.replace(".trim.aln.fasta", "").replace(".fasta", "")
        
        # 2. Get the specific model (default to WAG if missing from text file)
        base_model = og_models.get(og_name, "WAG")
        
        print(f"  Processing {og_name} (Model: {base_model})...")
        
        try:
            taxa_seqs = read_alignment(aln_path)

            # 3. Load IQ-TREE ML tree as BEAST2 starting tree (from Phase 1 -m MFP)
            starting_newick = read_iqtree_starting_tree(og_name)

            # 4. Generate two runs per orthogroup
            for run in ["run1", "run2"]:
                out_path = XML_DIR / f"{og_name}_{run}.xml"

                if out_path.exists():
                    print(f"    Skipping {og_name} {run} — XML already exists")
                    continue

                generate_xml(og_name, taxa_seqs, base_model, run, out_path, starting_newick)
                print(f"    Written: {out_path.name}  ({len(taxa_seqs)} taxa)")
                
        except Exception as e:
            print(f"  ERROR processing {og_name}: {e}")

    print()
    print(f"Done. XML files written to {XML_DIR}/")
    print()
    print("Next step — submit BEAST2 jobs:")
    print("  condor_submit scripts/beast_jobs.sub")


if __name__ == "__main__":
    main()