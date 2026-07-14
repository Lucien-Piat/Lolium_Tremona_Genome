#!/usr/bin/env python3
"""Assemble a composite figure: A = Circos, B = gene-proximal enrichment, C = composition donut."""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba

from plot_circos import build_plot

# Shared TE taxonomy and colors (kept consistent across A, B and C)
def category(c):
    if c == "LTR/Gypsy": return "Gypsy"
    if c == "LTR/Copia": return "Copia"
    if c.startswith("LTR/"): return "other LTR"
    if c.startswith("LINE"): return "LINE"
    if c.startswith("SINE"): return "SINE"
    if c.startswith("DNA/") or c.startswith("RC/"): return "DNA transposon"
    if c == "Unknown": return "Unknown"
    return "Other"

CATS = ["Gypsy", "Copia", "other LTR", "LINE", "DNA transposon", "SINE", "Unknown", "Other"]
CCOL = {"Gypsy": "#c0392b", "Copia": "#27ae60", "other LTR": "#16a085", "LINE": "#8e44ad",
        "DNA transposon": "#2980b9", "SINE": "#d2b4de", "Unknown": "#95a5a6", "Other": "#bdc3c7"}
ALPHA = 0.85
FACE = {c: to_rgba(CCOL[c], ALPHA) for c in CATS}
NEAR = 2000

CODING_COLOR = "#f1c40f"
UNANNOT_COLOR = "#ffffff"


def draw_enrichment(ax, dist_tsv):
    """Panel B: gene-proximal log2 enrichment by TE class (only plot kept from the old B)."""
    df = pd.read_csv(dist_tsv, sep="\t", names=["distance", "klass", "length"])
    df["cat"] = df.klass.map(category)

    allbp = df.groupby("cat")["length"].sum()
    nearbp = df[df.distance <= NEAR].groupby("cat")["length"].sum().reindex(allbp.index).fillna(0)
    enr = np.log2(((nearbp / nearbp.sum()) / (allbp / allbp.sum()))
                  .replace(0, np.nan)).dropna().sort_values()

    cols = [to_rgba(CCOL.get(c, "#999999"), ALPHA) for c in enr.index]
    ax.barh(range(len(enr)), enr.values, color=cols, edgecolor="black", linewidth=0.2)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(range(len(enr)))
    ax.set_yticklabels(enr.index, fontsize=9)
    ax.yaxis.tick_right()                       # ticks and labels on the right
    ax.yaxis.set_label_position("right")
    ax.set_xlabel("log2 enrichment within 2 kb of a gene")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def draw_composition_donut(ax, classtab, part_tsv):
    """Panel C: donut, bold names on short leaders outside, Gb size inside each wedge."""
    d = pd.read_csv(classtab, sep="\t")
    d.columns = ["klass", "frag", "bp", "pct"]
    d["cat"] = d.klass.map(category)

    p = pd.read_csv(part_tsv, sep="\t")
    part = dict(zip(p["component"], p["bp"]))
    genome = float(part["genome_total"])
    coding = float(part["coding_nonTE"])

    grouped = d.groupby("cat")["bp"].sum()
    te_total = float(d.bp.sum())

    segs = []
    for c in CATS:
        if c in grouped.index and grouped[c] > 0:
            segs.append((c, float(grouped[c]), CCOL[c]))
    segs.append(("Unannotated", max(genome - te_total - coding, 0.0), UNANNOT_COLOR))
    segs.append(("Coding DNA", coding, CODING_COLOR))

    sizes = [s[1] for s in segs]
    colors = [s[2] for s in segs]
    total_gb = sum(sizes) / 1e9

    wedges, _ = ax.pie(
        sizes, colors=colors, startangle=90, counterclock=False,
        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=2.0),
    )

    # ring midpoint radius (wedge spans 0.58 to 1.0 with width=0.42)
    r_in = 0.79
    arrow = dict(arrowstyle="-", color="#555555", lw=1.0)
    for wedge, (lab, bp, col) in zip(wedges, segs):
        gb = bp / 1e9
        ang = (wedge.theta2 - wedge.theta1) / 2.0 + wedge.theta1
        x = np.cos(np.deg2rad(ang))
        y = np.sin(np.deg2rad(ang))

        # bold name outside, pulled in close
        ha = "left" if x >= 0 else "right"
        conn = f"angle,angleA=0,angleB={ang}"
        ax.annotate(
            lab,
            xy=(x, y), xytext=(1.12 * np.sign(x), 1.05 * y),
            ha=ha, va="center", fontsize=13, fontweight="bold",
            arrowprops=dict(connectionstyle=conn, **arrow),
        )

        # Gb size inside the wedge
        if gb >= 0.05:
            txt_col = "white" if col not in (UNANNOT_COLOR, CODING_COLOR) else "#333"
            ax.text(r_in * x, r_in * y, f"{gb:.2f}", ha="center", va="center",
                    fontsize=10, fontweight="bold", color=txt_col)

    ax.set(aspect="equal")
    ax.text(0, 0, f"{total_gb:.2f}\nGb", ha="center", va="center",
            fontsize=13, fontweight="bold")

def draw_legend(ax):
    """Single shared legend below B and C."""
    items = [(c, FACE[c]) for c in CATS]
    items += [("Unannotated", to_rgba(UNANNOT_COLOR, 1.0)),
              ("Coding DNA", to_rgba(CODING_COLOR, 1.0))]
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=col, edgecolor="black", linewidth=0.3)
               for _, col in items]
    labels = [name for name, _ in items]
    ax.legend(handles, labels, loc="center", ncol=3, frameon=False,
              fontsize=9, handlelength=1.3, columnspacing=1.4, labelspacing=0.7)
    ax.axis("off")


def panel_letter(ax, text, x, y):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=20,
            fontweight="bold", va="bottom", ha="right")


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    fig = plt.figure(figsize=(22, 13))
    # Wider A, tighter gutter so the circos really takes the space
    outer = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.0], wspace=0.02)

    # A: circos on a polar axis in the big left column
    ax_circos = fig.add_subplot(outer[0, 0], projection="polar")

    # Right column: B (enrichment), C (donut), legend
    right = outer[0, 1].subgridspec(2, 1, height_ratios=[1.0, 1.3], hspace=0.4)
    axB = fig.add_subplot(right[0])
    axC = fig.add_subplot(right[1])

    print("[info] drawing panel A (circos)", file=sys.stderr)
    build_plot(args, ax=ax_circos)

    print("[info] drawing panel B (enrichment)", file=sys.stderr)
    draw_enrichment(axB, args.dist)

    print("[info] drawing panel C (composition donut)", file=sys.stderr)
    draw_composition_donut(axC, args.classtab, args.partition)


    panel_letter(ax_circos, "A.", x=0.02, y=0.98)
    panel_letter(axB, "B.", x=-0.18, y=1.02)
    panel_letter(axC, "C.", x=-0.05, y=1.02)

    fig.savefig(args.output, bbox_inches="tight", dpi=350)
    print(f"[info] wrote {args.output}", file=sys.stderr)


def parse_args():
    p = argparse.ArgumentParser()
    # Circos inputs (same names build_plot expects)
    p.add_argument("--genome", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--gff", required=True)
    p.add_argument("--busco", required=True)
    p.add_argument("--synteny", required=True)
    p.add_argument("--numt", required=True)
    p.add_argument("--nupt", required=True)
    p.add_argument("--mito-fasta", required=True, dest="mito_fasta")
    p.add_argument("--pltd-fasta", required=True, dest="pltd_fasta")
    p.add_argument("--mito-gb", required=True, dest="mito_gb")
    p.add_argument("--pltd-gb", required=True, dest="pltd_gb")
    p.add_argument("--te-gff", required=True, dest="te_gff")
    p.add_argument("--te-mapping", required=False, dest="te_mapping")
    # Panel B and C inputs
    p.add_argument("--dist", required=True, help="te_gene_distance.tsv for panel B")
    p.add_argument("--classtab", required=True, help="TE class_table for panel C")
    p.add_argument("--partition", required=True, help="genome_partition.tsv for panel C")
    # Composite output
    p.add_argument("--output", default="results/data_circo/circos_composite.pdf")
    return p.parse_args()


if __name__ == "__main__":
    main()