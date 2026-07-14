"""Build the Circos figure for Lolium multiflorum var Tremona."""
import argparse
import re
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib.cm as cm  # type: ignore
import numpy as np  # type: ignore
from pycirclize import Circos  # type: ignore

from tracks import (
    read_fai, organelle_lengths, find_gaps,
    gc_windows, organelle_gc, gene_density, busco_orthologs, organelle_features,
    load_self_synteny, load_organelle_links, load_te_mapping, te_density
)

ORG_SCALE = 150
GAP_MIN   = 10
WINDOW    = 1_000_000
ORG_WIN   = 2_000
ORG_SEP   = 5000

MIN_SYNTENY_GENES = 0
MIN_NUMT = 3000
MIN_NUPT = 2000

# Colors
CONTIG_GRAYS = ["#5a5a5a", "#a8a8a8"]
NUC_PALETTE  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf"]
MITO_COLOR   = "#8d6e63"
PLTD_COLOR   = "#7cb342"

TRACK_LABEL_SIZE = 10   
CHR_NAME_SIZE    = 17   
CHR_TICK_SIZE    = 11   

START_GAP = 15

BUSCO_COLORS = {
    "Complete":   "#2e7d32",
    "Duplicated": "#ef6c00",
}

# Synteny colors
SYNTENY_INTER = "#1976d2"
SYNTENY_INTRA = "#9c27b0"
NUMT_COLOR  = "#c2185b"
NUPT_COLOR  = "#558b2f"
ORG_FEAT_COLORS = {"CDS": "#37474f", "tRNA": "#0288d1", "rRNA": "#c62828"}

# Requested TE Colors
TE_COLORS = {
    "Gypsy": "#c0392b",
    "Copia": "#27ae60",
    "other LTR": "#16a085",
    "LINE": "#8e44ad",
    "SINE": "#d2b4de",
    "DNA transposon": "#2980b9",
    "Unknown": "#95a5a6",
    "Other": "#bdc3c7"
}
# Order to stack the Global TE Track neatly
TE_STACK_ORDER = ["Gypsy", "Copia", "other LTR", "LINE", "SINE", "DNA transposon", "Other", "Unknown"]

def sector_color(name):
    if name == "mito": return MITO_COLOR
    if name == "pltd": return PLTD_COLOR
    m = re.search(r"\d+", name)
    if m:
        return NUC_PALETTE[(int(m.group()) - 1) % len(NUC_PALETTE)]
    return "#9e9e9e"

def virtual_sector(contigs):
    offsets = OrderedDict()
    cum = 0
    for cid, clen in contigs.items():
        offsets[cid] = cum
        cum += clen + ORG_SEP
    return offsets, max(cum - ORG_SEP, 0)

def chromosome_segments(name, length, gaps):
    cg = sorted((gs, ge) for c, gs, ge in gaps if c == name)
    segs, cursor = [], 0
    for gs, ge in cg:
        if gs > cursor:
            segs.append((cursor, gs, False))
        segs.append((gs, ge, True))
        cursor = ge
    if cursor < length:
        segs.append((cursor, length, False))
    return segs

def draw_nuclear_ideogram(track, name, length, gaps):
    contig_idx = 0
    for s, e, is_gap in chromosome_segments(name, length, gaps):
        if is_gap:
            continue
        track.rect(s, e, fc=CONTIG_GRAYS[contig_idx % 2], ec="none")
        contig_idx += 1

def draw_organelle_ideogram(track, contigs, offsets):
    items = list(contigs.items())
    for i, (cid, clen) in enumerate(items):
        off = offsets[cid]
        track.rect(off * ORG_SCALE, (off + clen) * ORG_SCALE,
                   fc=CONTIG_GRAYS[i % 2], ec="none")

def build_plot(args, ax=None):
    print("[info] reading karyotype", file=sys.stderr)
    nuc = read_fai(args.fai)
    mito_contigs = organelle_lengths(args.mito_fasta)
    pltd_contigs = organelle_lengths(args.pltd_fasta)

    mito_offsets, mito_len = virtual_sector(mito_contigs)
    pltd_offsets, pltd_len = virtual_sector(pltd_contigs)

    sectors = dict(nuc)
    if mito_len > 0: sectors["mito"] = mito_len * ORG_SCALE
    if pltd_len > 0: sectors["pltd"] = pltd_len * ORG_SCALE

    print(f"[info] {len(nuc)} nuclear, {len(mito_contigs)} mito contigs, "
          f"{len(pltd_contigs)} pltd contigs", file=sys.stderr)

    print("[info] computing tracks", file=sys.stderr)
    gaps  = find_gaps(args.genome, min_gap=GAP_MIN)
    gc    = gc_windows(args.genome, window=WINDOW)

    all_nuc_gc = [v for _, _, _, v in gc]
    if all_nuc_gc:
        GC_VMIN = min(all_nuc_gc) - 1.0
        GC_VMAX = max(all_nuc_gc) + 1.0
    else:
        GC_VMIN, GC_VMAX = 35.0, 55.0
    GC_BASELINE = 45.0

    genes = gene_density(args.gff, nuc, window=WINDOW)
    busco = busco_orthologs(args.busco)
    mito_feats = organelle_features(args.mito_gb)
    pltd_feats = organelle_features(args.pltd_gb)
    mito_gc = organelle_gc(args.mito_fasta, mito_offsets, ORG_SCALE, ORG_WIN)
    pltd_gc = organelle_gc(args.pltd_fasta, pltd_offsets, ORG_SCALE, ORG_WIN)

    te_mapping_dict = load_te_mapping(args.te_mapping) if args.te_mapping else {}
    tes = te_density(args.te_gff, nuc, te_mapping_dict, window=WINDOW) if args.te_gff else {}

    synteny = load_self_synteny(args.synteny)
    numts = load_organelle_links(args.numt) if Path(args.numt).exists() else []
    nupts = load_organelle_links(args.nupt) if Path(args.nupt).exists() else []

    known_nuc, known_mito, known_pltd = set(nuc), set(mito_offsets), set(pltd_offsets)
    synteny = [l for l in synteny if l["q_chr"] in known_nuc and l["t_chr"] in known_nuc]
    numts   = [l for l in numts   if l["nuc_chr"] in known_nuc and l["org_chr"] in known_mito]
    nupts   = [l for l in nupts   if l["nuc_chr"] in known_nuc and l["org_chr"] in known_pltd]

    synteny = [l for l in synteny if l["n_genes"] >= MIN_SYNTENY_GENES]
    numts   = [l for l in numts   if (l["nuc_end"] - l["nuc_start"]) >= MIN_NUMT]
    nupts   = [l for l in nupts   if (l["nuc_end"] - l["nuc_start"]) >= MIN_NUPT]

    # Dynamically scale tracks
    global_te_total_max = 10
    if tes:
        total_max_vals = []
        for chrom, (x, widths, y_stack) in tes.items():
            if y_stack:
                total_cov = np.sum(list(y_stack.values()), axis=0)
                total_max_vals.append(np.max(total_cov))

        if total_max_vals: global_te_total_max = max(10, np.max(total_max_vals))

    print("[info] building Circos", file=sys.stderr)

    spaces = [3] * (len(sectors) - 1) + [START_GAP]
    circos = Circos(sectors, space=spaces)

    for sector in circos.sectors:
        name = sector.name
        is_mito = name == "mito"
        is_pltd = name == "pltd"
        is_org  = is_mito or is_pltd
        sec_len = sectors[name]

        # Explicit Labels (larger font: TRACK_LABEL_SIZE)
        if name == "chr1":
            lbl = dict(x=0, size=TRACK_LABEL_SIZE, weight="bold", ha="right")
            sector.text("Scaffolds          ", r=101.0, **lbl)
            sector.text("Contigs         ",     r=96.5,  **lbl)
            sector.text("TEs        ",          r=89.0,  **lbl)
            sector.text("Genes         ",       r=81.0,  **lbl)
            sector.text("GC%        ",          r=73.0,  **lbl)
            sector.text("BUSCO         ",       r=65.0,  **lbl)

        # Chromosome color band with LAST TICK ONLY
        color_tr = sector.add_track((100, 103))
        color_tr.axis(fc=sector_color(name), ec="black", lw=0.3)
        if not is_org:
            # Place exactly one tick at sec_len (Compatibility fix for older pycirclize)
            color_tr.xticks(
                [sec_len],
                tick_length=3,
                label_size=CHR_TICK_SIZE,
                label_orientation="vertical",
                labels=[f"{sec_len/1_000_000:.0f}M"]
            )

        # Ideogram
        ideo_tr = sector.add_track((94, 100))
        ideo_tr.axis(fc="white", ec="black", lw=0.3)
        if is_mito:
            draw_organelle_ideogram(ideo_tr, mito_contigs, mito_offsets)
        elif is_pltd:
            draw_organelle_ideogram(ideo_tr, pltd_contigs, pltd_offsets)
        else:
            draw_nuclear_ideogram(ideo_tr, name, sec_len, gaps)

        sector.text(name, r=112, size=CHR_NAME_SIZE, weight="bold")

        # Inner tracks mapping
        if is_org:
            # Invisible dummy track to anchor the organelle links deeply
            sector.add_track((62, 68))
        else:
            # 1. TOTAL TEs Track (Unified Stacked Bar, 86 - 92)
            te_total_tr = sector.add_track((86, 92))
            te_total_tr.axis(fc="white", ec="black", lw=0.2)
            if name in tes:
                x, widths, y_stack = tes[name]
                x = np.asarray(x, dtype=float)
                # anchor first/last points to the sector ends
                x_pad = np.concatenate(([0.0], x, [sec_len]))
                bottom = np.zeros(len(x_pad))
                for fam in TE_STACK_ORDER:
                    if fam in y_stack:
                        y = np.asarray(y_stack[fam], dtype=float)
                        y_pad = np.concatenate(([y[0]], y, [y[-1]]))
                        top = bottom + y_pad
                        te_total_tr.fill_between(
                            x_pad, top, y2=bottom,
                            vmin=0, vmax=global_te_total_max + 1,
                            fc=TE_COLORS[fam], ec="none", alpha=1.0,
                        )
                        bottom = top

            # Gene density
            gd_tr = sector.add_track((78, 84))
            gd_tr.axis(fc="white", ec="black", lw=0.2)
            vals = [v for c, _, _, v in genes if c == name]
            if vals:
                gd_tr.heatmap(vals, cmap="Greens", vmin=0, vmax=max(max(vals), 1))

            # GC content
            gc_tr = sector.add_track((70, 76))
            gc_tr.axis(fc="white", ec="black", lw=0.2)
            cgc = [(s, e, v) for c, s, e, v in gc if c == name]
            if cgc:
                x = np.array([(s + e) / 2 for s, e, _ in cgc])
                y = np.array([v for _, _, v in cgc])

                if len(x):
                    y_clipped = np.clip(y, GC_VMIN, GC_VMAX)
                    y_above   = np.maximum(y_clipped, GC_BASELINE)
                    y_below   = np.minimum(y_clipped, GC_BASELINE)
                    y_base    = np.full_like(y_clipped, GC_BASELINE)

                    gc_tr.fill_between(x, y_above, y2=y_base, vmin=GC_VMIN, vmax=GC_VMAX, fc="#4caf50", alpha=0.8)
                    gc_tr.fill_between(x, y_below, y2=y_base, vmin=GC_VMIN, vmax=GC_VMAX, fc="#e53935", alpha=0.8)
                    gc_tr.line(x, y_clipped, vmin=GC_VMIN, vmax=GC_VMAX, color="black", lw=0.3)

            # BUSCO
            bu_tr = sector.add_track((62, 68))
            bu_tr.axis(fc="white", ec="black", lw=0.2)
            cbu = [(pos, st) for c, pos, st in busco if c == name and st in BUSCO_COLORS]
            if cbu:
                bar_w = max(sec_len / 500, 100_000)
                for pos, status in cbu:
                    left  = max(0, pos - bar_w / 2)
                    right = min(sec_len, pos + bar_w / 2)
                    if right > left:
                        bu_tr.rect(left, right, fc=BUSCO_COLORS[status], ec="none", alpha=0.5)

    # Self-synteny links
    synteny.sort(key=lambda l: abs(l["q_start"] - l["q_end"]), reverse=True)
    for l in synteny:
        color = SYNTENY_INTRA if l["q_chr"] == l["t_chr"] else SYNTENY_INTER
        alpha = min(0.2 + l["n_genes"] / 50.0, 0.6)
        circos.link(
            (l["q_chr"], l["q_start"], l["q_end"]),
            (l["t_chr"], l["t_start"], l["t_end"]),
            color=color, alpha=alpha, lw=0.4,
        )

    numts.sort(key=lambda l: abs(l["nuc_start"] - l["nuc_end"]), reverse=True)
    for l in numts:
        off = mito_offsets[l["org_chr"]]
        circos.link(
            ("mito", (off + l["org_start"]) * ORG_SCALE, (off + l["org_end"]) * ORG_SCALE),
            (l["nuc_chr"], l["nuc_start"], l["nuc_end"]),
            color=NUMT_COLOR, alpha=0.6, lw=0.5,
        )

    nupts.sort(key=lambda l: abs(l["nuc_start"] - l["nuc_end"]), reverse=True)
    for l in nupts:
        off = pltd_offsets[l["org_chr"]]
        circos.link(
            ("pltd", (off + l["org_start"]) * ORG_SCALE, (off + l["org_end"]) * ORG_SCALE),
            (l["nuc_chr"], l["nuc_start"], l["nuc_end"]),
            color=NUPT_COLOR, alpha=0.6, lw=0.5,
        )

    # When an axis is supplied, draw onto it (composite figure) and return.
    if ax is not None:
        circos.plotfig(ax=ax)
        return circos

    fig = circos.plotfig(figsize=(14, 14), dpi=150)
    fig.savefig(args.output, bbox_inches="tight", dpi=450)
    print(f"[info] wrote {args.output}", file=sys.stderr)
    return fig

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--genome",      required=True)
    p.add_argument("--fai",         required=True)
    p.add_argument("--gff",         required=True)
    p.add_argument("--busco",       required=True)
    p.add_argument("--synteny",     required=True)
    p.add_argument("--numt",        required=True)
    p.add_argument("--nupt",        required=True)
    p.add_argument("--mito-fasta",  required=True, dest="mito_fasta")
    p.add_argument("--pltd-fasta",  required=True, dest="pltd_fasta")
    p.add_argument("--mito-gb",     required=True, dest="mito_gb")
    p.add_argument("--pltd-gb",     required=True, dest="pltd_gb")
    p.add_argument("--te-gff",      required=True, help="GFF3 file containing TE annotations")
    p.add_argument("--te-mapping",  required=False, help="TSV file mapping Motif ID to TE Class")
    p.add_argument("--output",      default="results/data_circo/circos.pdf")
    args = p.parse_args()
    build_plot(args)