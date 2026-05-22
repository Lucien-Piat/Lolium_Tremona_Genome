"""Build the Circos figure for Lolium multiflorum var Tremona."""
import argparse
import re
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib.cm as cm # type: ignore
import numpy as np # type: ignore
from matplotlib.patches import Patch # type: ignore
from pycirclize import Circos # type: ignore

from tracks import (
    read_fai, organelle_lengths, find_gaps,
    gc_windows, organelle_gc, gene_density, busco_orthologs, organelle_features,
    load_self_synteny, load_organelle_links,
)

# Visual parameters
ORG_SCALE = 200
GAP_MIN   = 10
WINDOW    = 1_000_000
ORG_WIN   = 2_000
ORG_SEP   = 5000

# GC range used across all sectors so nuclear and organelle scales are comparable
GC_VMIN, GC_VMAX, GC_BASELINE = 25.0, 65.0, 45.0

# Colors
CONTIG_GRAYS = ["#5a5a5a", "#a8a8a8"]
NUC_PALETTE  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf"]
MITO_COLOR   = "#8d6e63"
PLTD_COLOR   = "#7cb342"

BUSCO_COLORS = {
    "Complete":   "#2e7d32",  # green: complete single-copy
    "Duplicated": "#ef6c00",  # orange: complete duplicated
}

SYNTENY_FWD = "#1976d2"
SYNTENY_INV = "#fbc02d"
NUMT_COLOR  = "#c2185b"
NUPT_COLOR  = "#558b2f"

ORG_FEAT_COLORS = {"CDS": "#37474f", "tRNA": "#0288d1", "rRNA": "#c62828"}


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
            continue                            # no more red markers
        track.rect(s, e, fc=CONTIG_GRAYS[contig_idx % 2], ec="none")
        contig_idx += 1


def draw_organelle_ideogram(track, contigs, offsets):
    items = list(contigs.items())
    for i, (cid, clen) in enumerate(items):
        off = offsets[cid]
        track.rect(off * ORG_SCALE, (off + clen) * ORG_SCALE,
                   fc=CONTIG_GRAYS[i % 2], ec="none")


def build_plot(args):
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
    genes = gene_density(args.gff, nuc, window=WINDOW)
    busco = busco_orthologs(args.busco)
    mito_feats = organelle_features(args.mito_gb)
    pltd_feats = organelle_features(args.pltd_gb)
    mito_gc = organelle_gc(args.mito_fasta, mito_offsets, ORG_SCALE, ORG_WIN)
    pltd_gc = organelle_gc(args.pltd_fasta, pltd_offsets, ORG_SCALE, ORG_WIN)

    synteny = load_self_synteny(args.synteny)
    numts = load_organelle_links(args.numt) if Path(args.numt).exists() else []
    nupts = load_organelle_links(args.nupt) if Path(args.nupt).exists() else []

    known_nuc, known_mito, known_pltd = set(nuc), set(mito_offsets), set(pltd_offsets)
    synteny = [l for l in synteny if l["q_chr"] in known_nuc and l["t_chr"] in known_nuc]
    numts   = [l for l in numts   if l["nuc_chr"] in known_nuc and l["org_chr"] in known_mito]
    nupts   = [l for l in nupts   if l["nuc_chr"] in known_nuc and l["org_chr"] in known_pltd]
    print(f"[info] {len(synteny)} synteny, {len(numts)} NUMTs, {len(nupts)} NUPTs",
          file=sys.stderr)

    print("[info] building Circos", file=sys.stderr)
    circos = Circos(sectors, space=3)

    for sector in circos.sectors:
        name    = sector.name
        is_mito = name == "mito"
        is_pltd = name == "pltd"
        is_org  = is_mito or is_pltd
        sec_len = sectors[name]

        # Chromosome color band (outermost)
        color_tr = sector.add_track((100, 103))
        color_tr.axis(fc=sector_color(name), ec="black", lw=0.3)

        # Ideogram with alternating gray contigs (no gap markers)
        ideo_tr = sector.add_track((94, 100))
        ideo_tr.axis(fc="white", ec="black", lw=0.3)
        if is_mito:
            draw_organelle_ideogram(ideo_tr, mito_contigs, mito_offsets)
        elif is_pltd:
            draw_organelle_ideogram(ideo_tr, pltd_contigs, pltd_offsets)
        else:
            draw_nuclear_ideogram(ideo_tr, name, sec_len, gaps)

        sector.text(name, r=109, size=11, weight="bold")
        ideo_tr.xticks(by=10_000_000, tick_length=2, outer=True, label_formatter=lambda v: f"{v/1e6:.0f}M")
        # Read coverage placeholder
        reads_tr = sector.add_track((86, 92))
        reads_tr.axis(fc="#f5f5f5", ec="black", lw=0.2)

        # Repeats placeholder
        rep_tr = sector.add_track((78, 84))
        rep_tr.axis(fc="#f5f5f5", ec="black", lw=0.2)

        # Gene density (nuclear) or annotations (organelle)
        gd_tr = sector.add_track((70, 76))
        gd_tr.axis(fc="white", ec="black", lw=0.2)
        if is_org:
            feats = mito_feats if is_mito else pltd_feats
            offs  = mito_offsets if is_mito else pltd_offsets
            for cid, gname, fs, fe, strand, ftype in feats:
                if cid not in offs:
                    continue
                if fe - fs > 15000:
                    continue
                off = offs[cid]
                left  = max(0, (off + fs) * ORG_SCALE)
                right = min(sec_len, (off + fe) * ORG_SCALE)
                if right > left:
                    gd_tr.rect(left, right,
                               fc=ORG_FEAT_COLORS.get(ftype, "#888"), ec="none")
        else:
            vals = [v for c, _, _, v in genes if c == name]
            if vals:
                gd_tr.heatmap(vals, cmap="Greens",
                              vmin=0, vmax=max(max(vals), 1))

        # GC content (nuclear and organelles, shared scale)
        gc_tr = sector.add_track((62, 68))
        gc_tr.axis(fc="white", ec="black", lw=0.2)
        if is_org:
            data = mito_gc if is_mito else pltd_gc
            if data:
                x = np.array([(s + e) / 2 for s, e, _ in data])
                y = np.array([v for _, _, v in data])
        else:
            cgc = [(s, e, v) for c, s, e, v in gc if c == name]
            if cgc:
                x = np.array([(s + e) / 2 for s, e, _ in cgc])
                y = np.array([v for _, _, v in cgc])
            else:
                x, y = np.array([]), np.array([])
        if len(x):
            if is_org:
                vmin = y.min() - 2
                vmax = y.max() + 2
                vbase = y.mean()
            else:
                vmin, vmax, vbase = GC_VMIN, GC_VMAX, GC_BASELINE

            gc_tr.fill_between(x, y, y2=np.full_like(y, vbase),
                               vmin=vmin, vmax=vmax,
                               fc="#5c6bc0", alpha=0.4)
            gc_tr.line(x, y, vmin=vmin, vmax=vmax,
                       color="#1a237e", lw=0.5)

        # BUSCO: overlapping wider bars, only Complete + Duplicated
        bu_tr = sector.add_track((54, 60))
        bu_tr.axis(fc="white", ec="black", lw=0.2)
        if not is_org:
            cbu = [(pos, st) for c, pos, st in busco if c == name and st in BUSCO_COLORS]
            if cbu:
                bar_w = max(sec_len / 500, 100_000)
                for pos, status in cbu:
                    left  = max(0, pos - bar_w / 2)
                    right = min(sec_len, pos + bar_w / 2)
                    if right > left:
                        bu_tr.rect(left, right,
                                   fc=BUSCO_COLORS[status], ec="none", alpha=0.5)

    # Self-synteny links
    synteny.sort(key=lambda l: abs(l["q_start"] - l["q_end"]), reverse=True)
    for l in synteny:
        color = SYNTENY_FWD if l["strand"] == "+" else SYNTENY_INV
        alpha = min(0.2 + l["n_genes"] / 50.0, 0.6)
        circos.link(
            (l["q_chr"], l["q_start"], l["q_end"]),
            (l["t_chr"], l["t_start"], l["t_end"]),
            color=color, alpha=alpha, lw=0.4,
        )

    for l in numts:
        off = mito_offsets[l["org_chr"]]
        circos.link(
            ("mito", (off + l["org_start"]) * ORG_SCALE,
                     (off + l["org_end"])   * ORG_SCALE),
            (l["nuc_chr"], l["nuc_start"], l["nuc_end"]),
            color=NUMT_COLOR, alpha=0.6, lw=0.5,
        )
    for l in nupts:
        off = pltd_offsets[l["org_chr"]]
        circos.link(
            ("pltd", (off + l["org_start"]) * ORG_SCALE,
                     (off + l["org_end"])   * ORG_SCALE),
            (l["nuc_chr"], l["nuc_start"], l["nuc_end"]),
            color=NUPT_COLOR, alpha=0.6, lw=0.5,
        )

    fig = circos.plotfig(figsize=(14, 14), dpi=150)

    legend = [
        Patch(fc=CONTIG_GRAYS[0], label="Contig (odd)"),
        Patch(fc=CONTIG_GRAYS[1], label="Contig (even)"),
        Patch(fc="#bdbdbd",       label="Read coverage (placeholder)"),
        Patch(fc="#bdbdbd",       label="Repeats (placeholder)"),
        Patch(fc=cm.Greens(0.7),  label="Gene density"),
        Patch(fc="#5c6bc0",       label="GC content"),
        Patch(fc=BUSCO_COLORS["Complete"],   label="BUSCO Complete"),
        Patch(fc=BUSCO_COLORS["Duplicated"], label="BUSCO Duplicated"),
        Patch(fc=SYNTENY_FWD,     label="Self-synteny (forward)"),
        Patch(fc=SYNTENY_INV,     label="Self-synteny (inverted)"),
        Patch(fc=NUMT_COLOR,      label="NUMT (mito insertion)"),
        Patch(fc=NUPT_COLOR,      label="NUPT (plastid insertion)"),
    ]
    fig.legend(
        handles=legend, loc="lower right",
        bbox_to_anchor=(1.05, 0.0),
        fontsize=8, frameon=True,
        title=f"Lolium multiflorum var Tremona\n(organelles scaled {ORG_SCALE}x)",
        title_fontsize=9,
    )

    fig.savefig(args.output, bbox_inches="tight")
    if args.output.endswith(".pdf"):
        png = args.output[:-4] + ".png"
        fig.savefig(png, bbox_inches="tight", dpi=200)
        print(f"[info] wrote {png}", file=sys.stderr)
    print(f"[info] wrote {args.output}", file=sys.stderr)


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
    p.add_argument("--output",      default="results/data_circo/circos.pdf")
    args = p.parse_args()
    build_plot(args)