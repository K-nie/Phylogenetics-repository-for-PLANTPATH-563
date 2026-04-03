#!/usr/bin/env Rscript
# Author:Benjamin Narh-Madey
# run_dstatistic.R
# Calculates phylogenetic signal D-statistic for binary mannose trait
# across gene trees produced by RAxML-NG and/or BEAST2.
#
# Method: Fritz & Purvis (2010) D-statistic via caper::phylo.d()
# Same approach as Opulente et al. / Harrison et al. Y1000+ paper
# Reference: Chavez et al. 2026
#
# CSV FORMAT (tab-separated, three columns):
#   Species              — formal species name (for species tree, future use)
#   Species_name_in_OGs  — name used in FASTA/gene trees (used HERE for matching)
#   Mannose_Class        — binary trait: 1=growth, 0=no growth
#
# USAGE:
#   Rscript scripts/run_dstatistic.R \
#     --trait data/mannose_traits.csv \
#     --treedir results/ \
#     --output results/dstatistic_results.tsv
#
# OUTPUT:
#   TSV with one row per gene tree:
#   orthogroup, tree_type, n_tree_tips, n_csv_species, n_matched,
#   n_unmatched, pct_matched, D, p_random, p_brownian, interpretation
#
# Author: Benjamin Narh-Madey

suppressPackageStartupMessages({
    library(ape)
    library(caper)
    library(picante)
})

# ── Parse arguments ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    params <- list(
        trait   = NULL,
        treedir = NULL,
        output  = NULL,
        seed    = 42,
        perms   = 1000
    )
    i <- 1
    while (i <= length(args)) {
        if      (args[i] == "--trait")   { params$trait   <- args[i+1]; i <- i+2 }
        else if (args[i] == "--treedir") { params$treedir <- args[i+1]; i <- i+2 }
        else if (args[i] == "--output")  { params$output  <- args[i+1]; i <- i+2 }
        else if (args[i] == "--seed")    { params$seed    <- as.integer(args[i+1]); i <- i+2 }
        else if (args[i] == "--perms")   { params$perms   <- as.integer(args[i+1]); i <- i+2 }
        else { i <- i+1 }
    }
    return(params)
}

params <- parse_args(args)

# ── Set defaults ──────────────────────────────────────────────────────────────
BASE_DIR <- "/home/glbrc.org/narhmadey/2.Y1000_267_phylogenetics/2_Plant_Path_563"

if (is.null(params$trait))   params$trait   <- file.path(BASE_DIR, "data", "mannose_traits.csv")
if (is.null(params$treedir)) params$treedir <- file.path(BASE_DIR, "results")
if (is.null(params$output))  params$output  <- file.path(BASE_DIR, "results", "dstatistic_results.tsv")

# ── Column names matching your CSV format ─────────────────────────────────────
# Species              — formal species name (for species tree, future use)
# Species_name_in_OGs  — name used in FASTA/gene trees (used HERE for matching)
# Mannose_Class        — binary trait 1/0
OG_NAME_COL <- "Species_name_in_OGs"   # matches gene tree tip labels
TRAIT_COL   <- "Mannose_Class"          # binary mannose trait
SPECIES_COL <- "Species"                # formal name (reporting only)

cat("=======================================================\n")
cat("  Phylogenetic Signal D-statistic Analysis\n")
cat("  Method: Fritz & Purvis (2010) via caper::phylo.d()\n")
cat("  Reference: Chavez et al. 2026\n")
cat("  Author: Benjamin Narh-Madey\n")
cat("=======================================================\n")
cat(sprintf("Trait file:    %s\n", params$trait))
cat(sprintf("Tree dir:      %s\n", params$treedir))
cat(sprintf("Output:        %s\n", params$output))
cat(sprintf("Seed:          %d\n", params$seed))
cat(sprintf("Permutations:  %d\n", params$perms))
cat(sprintf("Name column:   %s (matched to gene tree tips)\n", OG_NAME_COL))
cat(sprintf("Trait column:  %s\n", TRAIT_COL))
cat("\n")

# ── Load trait data ───────────────────────────────────────────────────────────
if (!file.exists(params$trait)) {
    stop(sprintf("ERROR: Trait file not found: %s\n       Place mannose CSV in data/", params$trait))
}

# Read as tab-separated (three-column format)
trait_data <- read.table(params$trait,
                         header       = TRUE,
                         sep          = "\t",
                         stringsAsFactors = FALSE,
                         quote        = "",
                         comment.char = "")

# Validate required columns exist
required_cols <- c(OG_NAME_COL, TRAIT_COL, SPECIES_COL)
missing_cols  <- required_cols[!required_cols %in% colnames(trait_data)]
if (length(missing_cols) > 0) {
    cat(sprintf("ERROR: Missing columns: %s\n", paste(missing_cols, collapse=", ")))
    cat(sprintf("       Found columns:   %s\n", paste(colnames(trait_data), collapse=", ")))
    stop("Column mismatch — check CSV format")
}

# Remove rows with missing trait values
n_before   <- nrow(trait_data)
trait_data <- trait_data[!is.na(trait_data[[TRAIT_COL]]), ]
n_removed  <- n_before - nrow(trait_data)
if (n_removed > 0) {
    cat(sprintf("WARNING: Removed %d rows with missing %s values\n\n", n_removed, TRAIT_COL))
}

# Validate binary trait
invalid_vals <- unique(trait_data[[TRAIT_COL]])[!unique(trait_data[[TRAIT_COL]]) %in% c(0, 1)]
if (length(invalid_vals) > 0) {
    stop(sprintf("ERROR: %s contains non-binary values: %s\n       Only 0 and 1 are allowed",
                 TRAIT_COL, paste(invalid_vals, collapse=", ")))
}

cat("── Trait Data Summary ──────────────────────────────\n")
cat(sprintf("  Total species in CSV:    %d\n", nrow(trait_data)))
cat(sprintf("  Mannose_Class = 1:       %d (growth)\n",    sum(trait_data[[TRAIT_COL]] == 1)))
cat(sprintf("  Mannose_Class = 0:       %d (no growth)\n", sum(trait_data[[TRAIT_COL]] == 0)))
cat("\n")

# ── Find tree files ───────────────────────────────────────────────────────────
raxml_trees <- list.files(params$treedir,
                          pattern    = "\\.raxml\\.support$",
                          full.names = TRUE,
                          recursive  = FALSE)
beast_trees <- list.files(file.path(params$treedir, "beast"),
                          pattern    = "\\.beast\\.tree$",
                          full.names = TRUE,
                          recursive  = FALSE)

all_trees <- c(
    setNames(raxml_trees,
             paste0(sub(".*/(OG[^.]+).*", "\\1", raxml_trees), "_ML")),
    setNames(beast_trees,
             paste0(sub(".*/(OG[^.]+).*", "\\1", beast_trees), "_Bayesian"))
)

if (length(all_trees) == 0) {
    stop(sprintf(
        "ERROR: No tree files found\n  RAxML: %s/*.raxml.support\n  BEAST: %s/beast/*.beast.tree",
        params$treedir, params$treedir))
}

cat(sprintf("── Found %d tree file(s) ───────────────────────────\n", length(all_trees)))
for (nm in names(all_trees)) cat(sprintf("  %s\n", nm))
cat("\n")

# ── Helper: interpret D value ─────────────────────────────────────────────────
interpret_d <- function(D, p_random, p_brownian) {
    if (is.na(D))                return("FAILED")
    if (p_random > 0.05)         return("No significant phylogenetic signal (p_random > 0.05)")
    if (D < 0)                   return("Highly clustered — stronger than Brownian motion")
    if (D < 0.5)                 return("Clustered — strong phylogenetic signal")
    if (D < 1.0)                 return("Moderate phylogenetic signal")
    if (abs(D - 1.0) < 0.1)     return("Random — no phylogenetic signal")
    return("Overdispersed — convergent evolution")
}

# ── Helper: name matching diagnostic ─────────────────────────────────────────
check_name_match <- function(csv_names, tree_tips, og_name) {
    matched   <- csv_names[csv_names %in% tree_tips]
    unmatched <- csv_names[!csv_names %in% tree_tips]
    pct       <- round(100 * length(matched) / length(csv_names), 1)

    cat(sprintf("  Name matching for %s:\n", og_name))
    cat(sprintf("    CSV species:        %d\n", length(csv_names)))
    cat(sprintf("    Tree tips:          %d\n", length(tree_tips)))
    cat(sprintf("    Matched:            %d (%.1f%%)\n", length(matched), pct))
    cat(sprintf("    Unmatched in CSV:   %d\n", length(unmatched)))

    if (length(unmatched) > 0 && length(unmatched) <= 10) {
        cat("    Unmatched names:\n")
        for (nm in unmatched) cat(sprintf("      - %s\n", nm))
    } else if (length(unmatched) > 10) {
        cat(sprintf("    First 5 unmatched: %s ...\n",
                    paste(head(unmatched, 5), collapse=", ")))
    }
    if (pct < 50) {
        cat("    WARNING: <50% matched — check Species_name_in_OGs vs FASTA headers\n")
    }
    cat("\n")
    return(list(matched=matched, unmatched=unmatched, pct=pct))
}

# ── Run D-statistic for each tree ─────────────────────────────────────────────
set.seed(params$seed)
results <- list()

for (tree_name in names(all_trees)) {
    tree_file <- all_trees[[tree_name]]
    og        <- sub("_(ML|Bayesian)$", "", tree_name)
    tree_type <- sub(".*_(ML|Bayesian)$", "\\1", tree_name)

    cat("=======================================================\n")
    cat(sprintf("Processing: %s (%s tree)\n", og, tree_type))
    cat("=======================================================\n")

    # ── Load tree ─────────────────────────────────────────────────────────────
    tree <- tryCatch(
        read.tree(tree_file),
        error = function(e) {
            cat(sprintf("  ERROR reading tree: %s\n\n", e$message))
            return(NULL)
        }
    )
    if (is.null(tree)) {
        results[[tree_name]] <- data.frame(
            orthogroup=og, tree_type=tree_type, n_tree_tips=NA,
            n_csv_species=nrow(trait_data), n_matched=0,
            n_unmatched=NA, pct_matched=0, D=NA,
            p_random=NA, p_brownian=NA,
            interpretation="FAILED — could not read tree file",
            stringsAsFactors=FALSE)
        next
    }

    # ── Make tree binary and ultrametric ──────────────────────────────────────
    if (!is.binary(tree)) {
        cat("  Note: converting multifurcating tree to binary\n")
        tree <- multi2di(tree)
    }
    if (!is.ultrametric(tree)) {
        cat("  Note: making tree ultrametric with chronos()\n")
        tree <- tryCatch(
            chronos(tree, quiet=TRUE),
            error = function(e) {
                cat(sprintf("  WARNING: chronos() failed — trying correlated model\n"))
                chronos(tree, model="correlated", quiet=TRUE)
            }
        )
    }

    n_tree_tips <- length(tree$tip.label)

    # ── Name matching diagnostic ──────────────────────────────────────────────
    csv_names   <- trait_data[[OG_NAME_COL]]
    match_info  <- check_name_match(csv_names, tree$tip.label, og)
    n_csv       <- nrow(trait_data)
    n_matched   <- length(match_info$matched)
    n_unmatched <- length(match_info$unmatched)
    pct_matched <- match_info$pct

    if (n_matched < 10) {
        cat(sprintf("  SKIPPING — only %d matched taxa (minimum 10 required)\n\n", n_matched))
        results[[tree_name]] <- data.frame(
            orthogroup=og, tree_type=tree_type, n_tree_tips=n_tree_tips,
            n_csv_species=n_csv, n_matched=n_matched,
            n_unmatched=n_unmatched, pct_matched=pct_matched,
            D=NA, p_random=NA, p_brownian=NA,
            interpretation="SKIPPED — too few matched taxa (<10)",
            stringsAsFactors=FALSE)
        next
    }

    # ── Match trait data to tree ──────────────────────────────────────────────
    named_trait <- setNames(as.integer(trait_data[[TRAIT_COL]]),
                            trait_data[[OG_NAME_COL]])

    matched_obj <- tryCatch(
        match.phylo.data(tree, named_trait),
        error = function(e) {
            cat(sprintf("  ERROR in match.phylo.data(): %s\n\n", e$message))
            return(NULL)
        }
    )
    if (is.null(matched_obj)) {
        results[[tree_name]] <- data.frame(
            orthogroup=og, tree_type=tree_type, n_tree_tips=n_tree_tips,
            n_csv_species=n_csv, n_matched=n_matched,
            n_unmatched=n_unmatched, pct_matched=pct_matched,
            D=NA, p_random=NA, p_brownian=NA,
            interpretation="FAILED — match.phylo.data() error",
            stringsAsFactors=FALSE)
        next
    }

    trait_matched <- matched_obj$data
    tree_matched  <- matched_obj$phy

    # ── Check trait variation ─────────────────────────────────────────────────
    if (length(unique(trait_matched)) < 2) {
        cat(sprintf("  SKIPPING — no trait variation after matching\n\n"))
        results[[tree_name]] <- data.frame(
            orthogroup=og, tree_type=tree_type, n_tree_tips=n_tree_tips,
            n_csv_species=n_csv, n_matched=n_matched,
            n_unmatched=n_unmatched, pct_matched=pct_matched,
            D=NA, p_random=NA, p_brownian=NA,
            interpretation="SKIPPED — no trait variation in matched species",
            stringsAsFactors=FALSE)
        next
    }

    cat(sprintf("  Trait distribution in matched species:\n"))
    cat(sprintf("    Mannose=1 (growth):    %d\n", sum(trait_matched == 1)))
    cat(sprintf("    Mannose=0 (no growth): %d\n", sum(trait_matched == 0)))

    # ── Build comparative data object ─────────────────────────────────────────
    comp_df <- data.frame(
        og_name = names(trait_matched),
        mannose = as.integer(trait_matched),
        row.names = names(trait_matched),
        stringsAsFactors = FALSE
    )

    comp_data <- tryCatch(
        comparative.data(
            phy       = tree_matched,
            data      = comp_df,
            names.col = "og_name",
            vcv       = TRUE,
            warn.dropped = FALSE
        ),
        error = function(e) {
            cat(sprintf("  ERROR in comparative.data(): %s\n\n", e$message))
            return(NULL)
        }
    )
    if (is.null(comp_data)) {
        results[[tree_name]] <- data.frame(
            orthogroup=og, tree_type=tree_type, n_tree_tips=n_tree_tips,
            n_csv_species=n_csv, n_matched=n_matched,
            n_unmatched=n_unmatched, pct_matched=pct_matched,
            D=NA, p_random=NA, p_brownian=NA,
            interpretation="FAILED — comparative.data() error",
            stringsAsFactors=FALSE)
        next
    }

    # ── Calculate D-statistic ─────────────────────────────────────────────────
    cat(sprintf("  Running phylo.d() with %d permutations...\n", params$perms))

    d_result <- tryCatch(
        phylo.d(
            data   = comp_data,
            binvar = mannose,
            permut = params$perms
        ),
        error = function(e) {
            cat(sprintf("  ERROR in phylo.d(): %s\n\n", e$message))
            return(NULL)
        }
    )
    if (is.null(d_result)) {
        results[[tree_name]] <- data.frame(
            orthogroup=og, tree_type=tree_type, n_tree_tips=n_tree_tips,
            n_csv_species=n_csv, n_matched=n_matched,
            n_unmatched=n_unmatched, pct_matched=pct_matched,
            D=NA, p_random=NA, p_brownian=NA,
            interpretation="FAILED — phylo.d() error",
            stringsAsFactors=FALSE)
        next
    }

    # ── Extract and report results ────────────────────────────────────────────
    D_val      <- d_result$DEstimate
    p_random   <- d_result$Pval1
    p_brownian <- d_result$Pval0
    interp     <- interpret_d(D_val, p_random, p_brownian)

    cat(sprintf("\n  ── RESULTS ───────────────────────────────────\n"))
    cat(sprintf("  D-statistic:    %.4f\n", D_val))
    cat(sprintf("  p(random):      %.4f %s\n", p_random,
                ifelse(p_random   < 0.05, "  * significant", "")))
    cat(sprintf("  p(Brownian):    %.4f %s\n", p_brownian,
                ifelse(p_brownian < 0.05, "  * significant", "")))
    cat(sprintf("  Interpretation: %s\n\n", interp))

    results[[tree_name]] <- data.frame(
        orthogroup    = og,
        tree_type     = tree_type,
        n_tree_tips   = n_tree_tips,
        n_csv_species = n_csv,
        n_matched     = n_matched,
        n_unmatched   = n_unmatched,
        pct_matched   = pct_matched,
        D             = round(D_val, 4),
        p_random      = round(p_random, 4),
        p_brownian    = round(p_brownian, 4),
        interpretation = interp,
        stringsAsFactors = FALSE
    )
}

# ── Write results ─────────────────────────────────────────────────────────────
if (length(results) == 0) {
    stop("ERROR: No results produced — check tree files and trait data\n")
}

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL

dir.create(dirname(params$output), recursive=TRUE, showWarnings=FALSE)
write.table(results_df, params$output, sep="\t", quote=FALSE, row.names=FALSE)

# ── Final summary ─────────────────────────────────────────────────────────────
cat("=======================================================\n")
cat("  FINAL RESULTS SUMMARY\n")
cat("=======================================================\n")
print(results_df[, c("orthogroup", "tree_type", "n_matched",
                      "pct_matched", "D", "p_random",
                      "p_brownian", "interpretation")],
      row.names=FALSE)
cat("\n")
cat(sprintf("Full results written to: %s\n\n", params$output))
cat("── D-statistic interpretation ──────────────────────\n")
cat("  D <  0   : Highly clustered (stronger than Brownian)\n")
cat("  D =  0   : Brownian motion (strong phylogenetic signal)\n")
cat("  0 < D < 1: Moderate phylogenetic signal\n")
cat("  D =  1   : Random (no phylogenetic signal)\n")
cat("  D >  1   : Overdispersed (convergent evolution)\n\n")
cat("── p-value guide ───────────────────────────────────\n")
cat("  p(random)   < 0.05 → trait is non-randomly distributed\n")
cat("  p(Brownian) < 0.05 → trait deviates from Brownian motion\n\n")
cat("── Matching guide ──────────────────────────────────\n")
cat("  n_matched   = species in both CSV and tree\n")
cat("  pct_matched = % of CSV species found in tree\n")
cat("  Low pct_matched → check Species_name_in_OGs vs FASTA headers\n")
