"""
Per-block protein identity for syntenic paralogs.
"""
import argparse
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

HEADER_RE = re.compile(r"## Alignment\s+(\d+):.*?(\S+)&(\S+)\s+(plus|minus)")

INTER_COLOR = "#1976d2"
INTRA_COLOR = "#9c27b0"
HAPLO_COLOR = "#e53935"


def parse_collinearity(path):
    blocks, current = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("## Alignment"):
                m = HEADER_RE.search(line)
                if m:
                    if current:
                        blocks.append(current)
                    current = {
                        "id": int(m.group(1)),
                        "q_chr": m.group(2),
                        "t_chr": m.group(3),
                        "orient": m.group(4),
                        "pairs": [],
                    }
            elif current and line and not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    current["pairs"].append((parts[1].strip(), parts[2].strip()))
    if current:
        blocks.append(current)
    return blocks


def load_identity(blast_path):
    identity = {}
    with open(blast_path) as fh:
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) < 12:
                continue
            q, s, pid = p[0], p[1], float(p[2])
            key = (q, s)
            if pid > identity.get(key, -1):
                identity[key] = pid
    return identity


def block_mean_identity(block, identity_map):
    pidents = []
    for g1, g2 in block["pairs"]:
        v = identity_map.get((g1, g2)) or identity_map.get((g2, g1))
        if v is not None:
            pidents.append(v)
    return float(np.mean(pidents)) if pidents else None


def main(args):
    blocks = parse_collinearity(args.collinearity)
    identity = load_identity(args.blast)

    intra, inter = [], []
    for b in blocks:
        mi = block_mean_identity(b, identity)
        if mi is None:
            continue
        (intra if b["q_chr"] == b["t_chr"] else inter).append(mi)

    intra = np.array(intra); inter = np.array(inter)
    print(f"[info] intra n={len(intra)}, median={np.median(intra):.1f}%; "
          f"inter n={len(inter)}, median={np.median(inter):.1f}%", file=sys.stderr)

    haplo_lo, haplo_hi = 100.0 - args.heterozygosity, 100.0

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "legend.fontsize": 9,
    })

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Haplotype-duplication zone
    ax.axvspan(haplo_lo, haplo_hi, color=HAPLO_COLOR, alpha=0.15, zorder=0)

    bins = np.arange(50, 102, 2)
    centres = (bins[:-1] + bins[1:]) / 2
    
    bw = (bins[1] - bins[0]) * 0.85 
    
    inter_counts, _ = np.histogram(inter, bins=bins)
    intra_counts, _ = np.histogram(intra, bins=bins)

    ax.bar(centres, inter_counts, width=bw,
           color=INTER_COLOR, alpha=0.75, edgecolor="white", linewidth=0.5,
           label=f"Inter-chromosomal (n = {len(inter)})", zorder=2)
           
    ax.bar(centres, intra_counts, width=bw, bottom=inter_counts,
           color=INTRA_COLOR, alpha=0.75, edgecolor="white", linewidth=0.5,
           label=f"Intra-chromosomal (n = {len(intra)})", zorder=2)

    ax.set_xlabel("Mean protein identity per block (%)")
    ax.set_ylabel("Number of synteny blocks")
    ax.set_xlim(50, 100)
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98))

    plt.tight_layout()
    plt.savefig(args.output, bbox_inches="tight")
    png = args.output.rsplit(".", 1)[0] + ".png"
    plt.savefig(png, bbox_inches="tight", dpi=300)
    print(f"[info] wrote {args.output} and {png}", file=sys.stderr)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--collinearity",   required=True)
    p.add_argument("--blast",          required=True)
    p.add_argument("--heterozygosity", type=float, default=3.27)
    p.add_argument("--output",         default="results/data_circo/synteny_divergence.pdf")
    main(p.parse_args())