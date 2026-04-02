#!/usr/bin/env python3
import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt

def safe_mkdir(p):
    os.makedirs(p, exist_ok=True)

def pick_col(df, candidates):
    """Return first candidate column that exists (case-insensitive)."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None

def clean_species_name(pathlike: str) -> str:
    # e.g. yHMPu5000026055_diutina_siamensis_170912.final.cds -> diutina_siamensis
    base = os.path.basename(str(pathlike))
    base = re.sub(r"\.(final\.)?cds(\.gz)?$", "", base)
    base = re.sub(r"\.(fa|fna|fasta|ffn|faa|pep)$", "", base)
    # remove leading project/sample IDs like yHMPu..., yHAB..., etc.
    base = re.sub(r"^[A-Za-z0-9]+_", "", base)
    # remove trailing date-like tokens _170912 etc (6 digits) and _SPADES etc (keep if you want)
    base = re.sub(r"_[0-9]{6}$", "", base)
    return base

def main():
    if len(sys.argv) < 3:
        print("Usage: plot_seqkit_qc.py <seqkit_stats.tsv> <out_dir>")
        sys.exit(1)

    in_tsv = sys.argv[1]
    out_dir = sys.argv[2]
    safe_mkdir(out_dir)

    df = pd.read_csv(in_tsv, sep=r"\s+", engine="python", dtype=str)
    # Convert numeric columns where possible
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="ignore")

    file_col = pick_col(df, ["file", "filename", "path"])
    if file_col is None:
        # seqkit usually uses 'file'
        raise ValueError(f"Could not find a file column in: {list(df.columns)}")

    df["species"] = df[file_col].apply(clean_species_name)

    # Common seqkit columns (may vary slightly by version)
    num_seqs_col = pick_col(df, ["num_seqs", "numseqs", "seqs", "num_seqs:"])
    sum_len_col  = pick_col(df, ["sum_len", "sumlen", "sum_len(bp)", "sum_len:"])
    min_len_col  = pick_col(df, ["min_len", "minlen"])
    avg_len_col  = pick_col(df, ["avg_len", "avglen", "mean_len", "meanlen"])
    max_len_col  = pick_col(df, ["max_len", "maxlen"])
    gc_col       = pick_col(df, ["gc", "gc%","gc_percent","gc_content"])

    # Helper to make a histogram safely
    def hist_plot(series, title, xlabel, out_png, bins=50):
        s = pd.to_numeric(series, errors="coerce").dropna()
        if s.empty:
            return
        plt.figure()
        plt.hist(s, bins=bins)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(out_png, dpi=200)
        plt.close()

    # Helper: scatter
    def scatter_plot(x, y, title, xlabel, ylabel, out_png):
        x = pd.to_numeric(x, errors="coerce")
        y = pd.to_numeric(y, errors="coerce")
        m = x.notna() & y.notna()
        if m.sum() == 0:
            return
        plt.figure()
        plt.scatter(x[m], y[m], s=10)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(out_png, dpi=200)
        plt.close()

    # Plot 1: total CDS length per species (sum_len)
    if sum_len_col:
        hist_plot(df[sum_len_col],
                  "Total CDS length per species",
                  "Total CDS length (bp)",
                  os.path.join(out_dir, "hist_total_cds_length_sum_len.png"),
                  bins=60)

    # Plot 2: number of CDS sequences per species (num_seqs)
    if num_seqs_col:
        hist_plot(df[num_seqs_col],
                  "Number of CDS sequences per species",
                  "Number of CDS sequences",
                  os.path.join(out_dir, "hist_num_cds_sequences_num_seqs.png"),
                  bins=60)

    # Plot 3: average CDS length per species (avg_len)
    if avg_len_col:
        hist_plot(df[avg_len_col],
                  "Average CDS length per species",
                  "Average CDS length (bp)",
                  os.path.join(out_dir, "hist_avg_cds_length_avg_len.png"),
                  bins=60)

    # Plot 4: min/max CDS length distributions (if present)
    if min_len_col:
        hist_plot(df[min_len_col],
                  "Minimum CDS length per species",
                  "Minimum CDS length (bp)",
                  os.path.join(out_dir, "hist_min_cds_length_min_len.png"),
                  bins=60)
    if max_len_col:
        hist_plot(df[max_len_col],
                  "Maximum CDS length per species",
                  "Maximum CDS length (bp)",
                  os.path.join(out_dir, "hist_max_cds_length_max_len.png"),
                  bins=60)

    # Plot 5: Scatter num_seqs vs sum_len (gene count vs total CDS bp)
    if num_seqs_col and sum_len_col:
        scatter_plot(df[num_seqs_col], df[sum_len_col],
                     "Num CDS sequences vs total CDS length",
                     "Num CDS sequences",
                     "Total CDS length (bp)",
                     os.path.join(out_dir, "scatter_num_seqs_vs_sum_len.png"))

    # Plot 6: GC% distribution (if present)
    if gc_col:
        hist_plot(df[gc_col],
                  "GC% distribution across species",
                  "GC%",
                  os.path.join(out_dir, "hist_gc_percent.png"),
                  bins=60)

    # Write an outlier report (top/bottom 10 by sum_len and num_seqs)
    report_path = os.path.join(out_dir, "qc_outliers_top_bottom.txt")
    with open(report_path, "w") as f:
        f.write("QC outlier summary\n\n")
        if sum_len_col:
            f.write(f"Top 10 by {sum_len_col}:\n")
            f.write(df.sort_values(sum_len_col, ascending=False)[["species", sum_len_col]].head(10).to_string(index=False))
            f.write("\n\nBottom 10 by {0}:\n".format(sum_len_col))
            f.write(df.sort_values(sum_len_col, ascending=True)[["species", sum_len_col]].head(10).to_string(index=False))
            f.write("\n\n")
        if num_seqs_col:
            f.write(f"Top 10 by {num_seqs_col}:\n")
            f.write(df.sort_values(num_seqs_col, ascending=False)[["species", num_seqs_col]].head(10).to_string(index=False))
            f.write("\n\nBottom 10 by {0}:\n".format(num_seqs_col))
            f.write(df.sort_values(num_seqs_col, ascending=True)[["species", num_seqs_col]].head(10).to_string(index=False))
            f.write("\n")

    print(f"[DONE] Plots written to: {out_dir}")
    print(f"[DONE] Outlier report: {report_path}")

if __name__ == "__main__":
    main()
