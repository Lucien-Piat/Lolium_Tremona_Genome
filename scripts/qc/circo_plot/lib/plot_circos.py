"""Build the Circos figure for Lolium multiflorum var Tremona."""
import argparse
import re
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib.cm as cm  # type: ignore
import numpy as np  # type: ignore
from matplotlib.patches import Patch  # type: ignore
from pycirclize import Circos  # type: ignore

from tracks import (
    read_fai, organelle_lengths, find_gaps,
    gc_windows, organelle_gc, gene_density, busco_orthologs, organelle_features,
    load_self_synteny, load_organelle_links, load_and_average_coverage
)

ORG_SCALE = 200
GAP_MIN   = 10
WINDOW    = 1_000_000
ORG_WIN   = 2_000
ORG_SEP   = 5000

MIN_SYNTENY_GENES = 10
MIN_NUMT = 3000
MIN_NUPT = 2000
COV_PCTILE        = 1

# Colors
CONTIG_GRAYS = ["#5a5a5a", "#a8a8a8"]
NUC_PALETTE  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf"]
MITO_COLOR   = "#8d6e63"
PLTD_COLOR   = "#7cb342"

BUSCO_COLORS = {
    "Complete":   "#2e7d32",
    "Duplicated": "#ef6c00",
}

# Synteny colors
SYNTENY_FWD   = "#1976d2"   # full variant: forward
SYNTENY_INV   = "#fbc02d"   # full variant: inverted
SYNTENY_INTER = "#1976d2"   # clean variant: between chromosomes
SYNTENY_INTRA = "#9c27b0"   # clean variant: within one chromosome

NUMT_COLOR  = "#c2185b"
NUPT_COLOR  = "#558b2f"

# Coverage anomaly markers
COV_MARK_LOW  = "#d32f2f"
COV_MARK_HIGH = "#fbc02d"

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
            continue
        track.rect(s, e, fc=CONTIG_GRAYS[contig_idx % 2], ec="none")
        contig_idx += 1


def draw_organelle_ideogram(track, contigs, offsets):
    items = list(contigs.items())
    for i, (cid, clen) in enumerate(items):
        off = offsets[cid]
        track.rect(off * ORG_SCALE, (off + clen) * ORG_SCALE,
                   fc=CONTIG_GRAYS[i % 2], ec="none")


def build_plot(args):
    is_clean = args.variant == "clean"
    print(f"[info] variant: {args.variant}", file=sys.stderr)

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
    print(f"[info] Auto GC scale (Nuclear): min={GC_VMIN:.2f}%, max={GC_VMAX:.2f}%", file=sys.stderr)

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

    synteny = [l for l in synteny if l["n_genes"] >= MIN_SYNTENY_GENES]
    numts   = [l for l in numts   if (l["nuc_end"] - l["nuc_start"]) >= MIN_NUMT]
    nupts   = [l for l in nupts   if (l["nuc_end"] - l["nuc_start"]) >= MIN_NUPT]
    print(f"[info] {len(synteny)} synteny (>= {MIN_SYNTENY_GENES} genes), "
          f"{len(numts)} NUMTs (>= {MIN_NUMT} bp), "
          f"{len(nupts)} NUPTs (>= {MIN_NUPT} bp)", file=sys.stderr)

    coverage_data = load_and_average_coverage(args.coverage)
    global_cov_max  = 100
    global_cov_min  = 0
    global_cov_base = 50
    cov_low_thr  = None
    cov_high_thr = None
    if coverage_data:
        all_cov_vals = [v for _, _, _, v in coverage_data]
        global_cov_max  = np.percentile(all_cov_vals, 98)
        global_cov_min  = max(0, np.percentile(all_cov_vals, 2))
        global_cov_base = np.median(all_cov_vals)
        cov_low_thr  = np.percentile(all_cov_vals, COV_PCTILE)
        cov_high_thr = np.percentile(all_cov_vals, 100 - COV_PCTILE)
        print(f"[info] Coverage anomaly thresholds: low<{cov_low_thr:.0f}, "
              f"high>{cov_high_thr:.0f}", file=sys.stderr)
    print(f"[info] Auto Coverage scale: min={global_cov_min:.0f}, max={global_cov_max:.0f}, "
          f"base={global_cov_base:.0f}", file=sys.stderr)

    print("[info] building Circos", file=sys.stderr)
    spaces = [3] * (len(sectors) - 1) + [5]
    circos = Circos(sectors, space=spaces)

    for sector in circos.sectors:
        name = sector.name

        if name == "chr1":
            sector.text("a.      ", r=97, x=0, size=11, weight="bold", ha="right")
            sector.text("b.      ", r=89, x=0, size=11, weight="bold", ha="right")
            sector.text("c.      ", r=81, x=0, size=11, weight="bold", ha="right")
            sector.text("d.      ", r=73, x=0, size=11, weight="bold", ha="right")
            sector.text("e.      ", r=65, x=0, size=11, weight="bold", ha="right")
            sector.text("f.      ", r=57, x=0, size=11, weight="bold", ha="right")

        is_mito = name == "mito"
        is_pltd = name == "pltd"
        is_org  = is_mito or is_pltd
        sec_len = sectors[name]

        # Chromosome color band
        color_tr = sector.add_track((100, 103))
        color_tr.axis(fc=sector_color(name), ec="black", lw=0.3)

        # Ideogram
        ideo_tr = sector.add_track((94, 100))
        ideo_tr.axis(fc="white", ec="black", lw=0.3)
        if is_mito:
            draw_organelle_ideogram(ideo_tr, mito_contigs, mito_offsets)
        elif is_pltd:
            draw_organelle_ideogram(ideo_tr, pltd_contigs, pltd_offsets)
        else:
            draw_nuclear_ideogram(ideo_tr, name, sec_len, gaps)

        sector.text(name, r=109, size=11, weight="bold")

        # Read coverage
        reads_tr = sector.add_track((86, 92))
        reads_tr.axis(fc="white", ec="black", lw=0.2)

        if not is_org and coverage_data:
            cov_sector = [(s, e, v) for c, s, e, v in coverage_data if c == name]
            if cov_sector:
                cov_sector.sort(key=lambda x: x[0])
                x = np.array([(s + e) / 2 for s, e, _ in cov_sector])
                y = np.array([v for _, _, v in cov_sector])

                y_clipped = np.clip(y, global_cov_min, global_cov_max)

                window_size = 20
                y_smoothed = np.convolve(y_clipped, np.ones(window_size)/window_size, mode='same')
                y_smoothed = np.clip(y_smoothed, global_cov_min, global_cov_max)
                y_base  = np.full_like(y_smoothed, global_cov_base)
                y_above = np.maximum(y_smoothed, global_cov_base)
                y_below = np.minimum(y_smoothed, global_cov_base)

                reads_tr.fill_between(x, y_above, y2=y_base,
                                      vmin=global_cov_min, vmax=global_cov_max,
                                      fc="#5e35b1", alpha=0.85)
                reads_tr.fill_between(x, y_below, y2=y_base,
                                      vmin=global_cov_min, vmax=global_cov_max,
                                      fc="#b39ddb", alpha=0.6)
                reads_tr.line(x, y_smoothed, vmin=global_cov_min, vmax=global_cov_max,
                              color="#311b92", lw=0.6)

                # Anomaly markers along top edge of coverage track
                if cov_low_thr is not None:
                    marker_w = max(sec_len / 1500, 30_000)
                    for s, e, v in cov_sector:
                        mid = (s + e) / 2
                        left  = max(0, mid - marker_w / 2)
                        right = min(sec_len, mid + marker_w / 2)
                        if right <= left:
                            continue
                        if v < cov_low_thr:
                            reads_tr.rect(left, right, r_lim=(91, 92),
                                          fc=COV_MARK_LOW, ec="none", alpha=0.95)
                        elif v > cov_high_thr:
                            reads_tr.rect(left, right, r_lim=(91, 92),
                                          fc=COV_MARK_HIGH, ec="none", alpha=0.95)

        # Repeats placeholder
        rep_tr = sector.add_track((78, 84))
        rep_tr.axis(fc="#f5f5f5", ec="black", lw=0.2)

        # Gene density / organelle annotations
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

        # GC content
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
            y_clipped = np.clip(y, GC_VMIN, GC_VMAX)
            y_above   = np.maximum(y_clipped, GC_BASELINE)
            y_below   = np.minimum(y_clipped, GC_BASELINE)
            y_base    = np.full_like(y_clipped, GC_BASELINE)

            gc_tr.fill_between(x, y_above, y2=y_base,
                               vmin=GC_VMIN, vmax=GC_VMAX,
                               fc="#4caf50", alpha=0.8)
            gc_tr.fill_between(x, y_below, y2=y_base,
                               vmin=GC_VMIN, vmax=GC_VMAX,
                               fc="#e53935", alpha=0.8)
            gc_tr.line(x, y_clipped, vmin=GC_VMIN, vmax=GC_VMAX,
                       color="black", lw=0.3)

        # BUSCO
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

    # Self-synteny links: only the color scheme differs between variants
    synteny.sort(key=lambda l: abs(l["q_start"] - l["q_end"]), reverse=True)
    for l in synteny:
        if is_clean:
            color = SYNTENY_INTRA if l["q_chr"] == l["t_chr"] else SYNTENY_INTER
        else:
            color = SYNTENY_FWD if l["strand"] == "+" else SYNTENY_INV
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
            ("mito", (off + l["org_start"]) * ORG_SCALE,
                     (off + l["org_end"])   * ORG_SCALE),
            (l["nuc_chr"], l["nuc_start"], l["nuc_end"]),
            color=NUMT_COLOR, alpha=0.6, lw=0.5,
        )

    nupts.sort(key=lambda l: abs(l["nuc_start"] - l["nuc_end"]), reverse=True)
    for l in nupts:
        off = pltd_offsets[l["org_chr"]]
        circos.link(
            ("pltd", (off + l["org_start"]) * ORG_SCALE,
                     (off + l["org_end"])   * ORG_SCALE),
            (l["nuc_chr"], l["nuc_start"], l["nuc_end"]),
            color=NUPT_COLOR, alpha=0.6, lw=0.5,
        )

    fig = circos.plotfig(figsize=(14, 14), dpi=150)

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
    p.add_argument("--coverage",    nargs="*", default=[], help="Coverage files (multiple allowed)")
    p.add_argument("--variant",     choices=["full", "clean"], default="full",
                   help="'full' (fwd/inv colors) or 'clean' (intra/inter colors)")
    p.add_argument("--output",      default="results/data_circo/circos.pdf")
    args = p.parse_args()
    build_plot(args)