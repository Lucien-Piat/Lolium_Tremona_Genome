#!/usr/bin/env python3
"""
  1. Read an assembly list:  fasta_path <TAB> busco_path <TAB> label
  2. Parse each FASTA, find N-runs, emit one row per contig.
  3. Parse each BUSCO full_table.tsv, map original sequence names to standard chroms.
  4. Plot one figure per chromosome with alternating grey contig blocks + BUSCO hits on the bottom third.

Usage:
  python pangenome_busco_plot.py assemblies_with_busco.tsv -n 7 -t 4
"""

import argparse
import gzip
import re
import sys
import concurrent.futures
from pathlib import Path

import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import matplotlib.patches as mpatches # type: ignore
from matplotlib.ticker import FuncFormatter, MultipleLocator # type: ignore

def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def iter_fasta(fh):
    name, chunks = None, []
    for line in fh:
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)

def find_contigs(seq, min_gap=10):
    gap_re = re.compile(f"N{{{min_gap},}}", re.IGNORECASE)
    out, last = [], 0
    for m in gap_re.finditer(seq):
        if m.start() > last:
            out.append((last + 1, m.start()))
        last = m.end()
    if last < len(seq):
        out.append((last + 1, len(seq)))
    if not out and len(seq) > 0:
        out = [(1, len(seq))]
    return out

def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]

def parse_list(path):
    entries = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"[warning] Skipping malformed line: {line}", file=sys.stderr)
                continue
            assembly = parts[0]
            busco = parts[1]
            label = parts[2] if len(parts) >= 3 else Path(assembly).stem
            entries.append((assembly, busco, label))
    return entries

def _process_single_assembly(path, label, n_chrom, min_gap, chrom_pattern):
    print(f"[info] reading FASTA for {label} <- {path}", file=sys.stderr)
    seqs = []
    with opener(path) as fh:
        for name, seq in iter_fasta(fh):
            seqs.append((name, seq.upper()))

    if chrom_pattern:
        pat = re.compile(chrom_pattern)
        picked = [(n, s) for n, s in seqs if pat.search(n)]
    else:
        seqs.sort(key=lambda x: -len(x[1]))
        picked = seqs[:n_chrom]

    picked.sort(key=lambda x: natural_key(x[0]))

    rows = []
    for i, (name, seq) in enumerate(picked, 1):
        canonical = f"chr{i:02d}"
        for j, (s, e) in enumerate(find_contigs(seq, min_gap=min_gap), 1):
            rows.append({
                "genome": label,
                "original_chrom": name,
                "chrom": canonical,
                "contig_idx": j,
                "start": s,
                "end": e,
                "length": e - s + 1,
            })
    return rows

def extract_contigs(entries, n_chrom, min_gap=10, chrom_pattern=None, threads=1):
    rows = []
    if threads <= 1:
        for path, _, label in entries:
            rows.extend(_process_single_assembly(path, label, n_chrom, min_gap, chrom_pattern))
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [
                executor.submit(_process_single_assembly, path, label, n_chrom, min_gap, chrom_pattern)
                for path, _, label in entries
            ]
            for future in futures:
                rows.extend(future.result())
    return pd.DataFrame(rows)


# ============================ BUSCO parsing ==================================

def load_busco(busco_path, label, chrom_map):
    print(f"[info] reading BUSCO for {label} <- {busco_path}", file=sys.stderr)
    
    with open(busco_path) as f:
        lines = f.readlines()
        
    header = None
    data = []
    for line in lines:
        if line.startswith("# Busco id") or line.startswith("# Busco ID"):
            header = line.strip().lstrip("#").split("\t")
        elif not line.startswith("#") and line.strip():
            data.append(line.strip().split("\t"))
            
    if not header:
        header = ["BuscoId", "Status", "Sequence", "Start", "End"]
        
    df = pd.DataFrame(data, columns=header[:len(data[0])] if data else header)
    df.columns = [c.strip() for c in df.columns]
    
    rename_dict = {"Busco id": "BuscoId", "Gene Start": "Start", "Gene End": "End"}
    df = df.rename(columns={k: v for k, v in rename_dict.items() if k in df.columns})
    
    df["Status"] = df["Status"].replace({"Complete": "Single"})
    df["genome"] = label
    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    
    df["chrom"] = df["Sequence"].map(chrom_map)
    df = df.dropna(subset=["chrom", "Start"])
    
    return df


# VISUAL CONFIGURATION
DARK = "#959595"
LIGHT = "#dcdcdc"
BUSCO_COLORS = {"Single": "#4CAF50", "Duplicated": "#FF9800", "Fragmented": "#F44336"}

BUSCO_LINEWIDTH = 0.9

A4_W = 10
A4_H = 11.69

def plot_a4(df_contigs, df_busco, genome_order, chroms, tick_step_mb=10, bar_height=0.7):
    n_chr = len(chroms)
    n_gen = len(genome_order)

    fig, axes = plt.subplots(n_chr, 1, figsize=(A4_W, A4_H), gridspec_kw={"hspace": 0.55})
    axes = [axes] if n_chr == 1 else axes 

    y_pos = {g: n_gen - 1 - i for i, g in enumerate(genome_order)}

    for ax, chrom in zip(axes, chroms):
        df_chr = df_contigs[df_contigs["chrom"] == chrom]
        
        if df_chr.empty:
            continue
            
        max_len = df_chr.groupby("genome")["end"].max().max()

        for g in genome_order:
            y = y_pos[g]
            
            # 1. Contig blocks (Grey)
            sub_contig = df_chr[df_chr["genome"] == g].sort_values("start")
            if sub_contig.empty:
                continue
                
            for _, r in sub_contig.iterrows():
                color = DARK if r["contig_idx"] % 2 == 1 else LIGHT
                ax.add_patch(mpatches.Rectangle(
                    (r["start"], y - bar_height / 2),
                    r["end"] - r["start"] + 1, bar_height,
                    facecolor=color, edgecolor="none", linewidth=0))
            
            total_end = sub_contig["end"].max()
            
            ax.add_patch(mpatches.Rectangle(
                (1, y - bar_height / 2), total_end - 1, bar_height,
                facecolor="none", edgecolor="black", linewidth=0.4))

            # 2. BUSCO vlines
            sub_busco = df_busco[(df_busco["chrom"] == chrom) & (df_busco["genome"] == g)]
            
            # Coordinates of the vlines
            busco_ymin = y - bar_height / 2
            busco_ymax = busco_ymin + (bar_height / 3)
            
            for status, col in BUSCO_COLORS.items():
                pts = sub_busco.loc[sub_busco["Status"] == status, "Start"]
                if len(pts):
                    ax.vlines(pts, busco_ymin, busco_ymax,
                              colors=col, linewidth=BUSCO_LINEWIDTH, alpha=0.8)

        ax.set_xlim(0, max_len * 1.02)
        ax.set_ylim(-0.5, n_gen - 0.5)
        ax.set_yticks(list(y_pos.values()))
        ax.set_yticklabels(list(y_pos.keys()), fontsize=9)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, p: f"{x/1e6:.0f}"))
        ax.xaxis.set_major_locator(MultipleLocator(tick_step_mb * 1e6))
        ax.tick_params(axis="x", labelsize=5)
        ax.tick_params(axis="y", length=0)
        ax.set_title(chrom, fontsize=7, loc="left", pad=2)
        for s in ["top", "right", "left"]:
            ax.spines[s].set_visible(False)

    axes[-1].set_xlabel("Position (Mb)", fontsize=13)
    
    fig.subplots_adjust(left=0.08, right=0.96, top=0.97, bottom=0.04)
    return fig

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("assembly_list", help="TSV: fasta_path <TAB> busco_path <TAB> label")
    ap.add_argument("-n", "--n-chromosomes", type=int, required=True, help="Number of chromosomes per assembly")
    ap.add_argument("-t", "--threads", type=int, default=1, help="Number of CPU threads for parsing")
    ap.add_argument("--min-gap", type=int, default=10, help="Minimum N-run length to treat as scaffold gap")
    ap.add_argument("--chrom-pattern", default=None, help="Regex to select chromosome sequence names")
    ap.add_argument("--contigs-tsv", default="contigs_all.tsv")
    ap.add_argument("--output-pdf", default="pangenome_busco_chromosomes.pdf")
    ap.add_argument("--tick-step-mb", type=float, default=10.0)
    args = ap.parse_args()

    entries = parse_list(args.assembly_list)
    print(f"[info] {len(entries)} assemblies to process", file=sys.stderr)

    # 1. Parse Contigs
    df_contigs = extract_contigs(entries, args.n_chromosomes, min_gap=args.min_gap,
                                 chrom_pattern=args.chrom_pattern, threads=args.threads)
    df_contigs.to_csv(args.contigs_tsv, sep="\t", index=False)
    
    chrom_map = {}
    for g in df_contigs["genome"].unique():
        sub = df_contigs[df_contigs["genome"] == g]
        chrom_map[g] = dict(zip(sub["original_chrom"], sub["chrom"]))

    # 2. Parse BUSCO
    busco_dfs = []
    for _, busco_path, label in entries:
        b_df = load_busco(busco_path, label, chrom_map.get(label, {}))
        busco_dfs.append(b_df)
    
    df_busco = pd.concat(busco_dfs, ignore_index=True) if busco_dfs else pd.DataFrame()

    genome_order = [label for _, _, label in entries]
    chroms = sorted(df_contigs["chrom"].unique(), key=natural_key)

    # 3. Draw Plot
    fig = plot_a4(df_contigs, df_busco, genome_order, chroms, tick_step_mb=args.tick_step_mb)
    
    fig.savefig(args.output_pdf, format="pdf")
    plt.close(fig)
    print(f"[info] wrote {args.output_pdf}", file=sys.stderr)

if __name__ == "__main__":
    main()