"""Build the Circos figure for Lolium multiflorum var Tremona."""
import argparse
import sys
from pathlib import Path

import matplotlib as mpl # type: ignore
import matplotlib.cm as cm # type: ignore
import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore
from matplotlib.patches import Patch # type: ignore
from pycirclize import Circos # type: ignore

from tracks import (
    read_fai, organelle_lengths, find_gaps,
    gc_windows, gene_density, busco_orthologs, organelle_features,
    load_self_synteny, load_organelle_links,
)

# Visual parameters
ORG_SCALE = 200            # scale organelles up so they're visible
GAP_MIN  = 1000
WINDOW   = 1_000_000

BUSCO_COLORS = {
    "Complete":   "#2e7d32",
    "Duplicated": "#1565c0",
    "Fragmented": "#ef6c00",
}
BUSCO_Y = {"Complete": 0.8, "Duplicated": 0.5, "Fragmented": 0.2}

SYNTENY_FWD = "#1976d2"
SYNTENY_INV = "#d32f2f"
NUMT_COLOR  = "#c2185b"
NUPT_COLOR  = "#558b2f"
GAP_COLOR   = "#b71c1c"

ORG_FEAT_COLORS = {"CDS": "#37474f", "tRNA": "#0288d1", "rRNA": "#c62828"}


def chromosome_palette(names):
    cmap = cm.get_cmap("Set3", max(len(names), 3))
    return {n: cmap(i) for i, n in enumerate(names)}


def organelle_color(name):
    low = name.lower()
    if any(k in low for k in ("mito", "mt")):
        return "#8d6e63"
    if any(k in low for k in ("pltd", "chloro", "plast", "cp")):
        return "#7cb342"
    return "#9e9e9e"


def build_plot(args):
    nuc = read_fai(args.fai)
    org = organelle_lengths(args.mito_fasta, args.pltd_fasta)

    sectors = dict(nuc)
    for name, length in org.items():
        sectors[name] = length * ORG_SCALE

    print(f"[info] {len(nuc)} nuclear + {len(org)} organelle sectors",
          file=sys.stderr)

    gaps  = find_gaps(args.genome, min_gap=GAP_MIN)
    gc    = gc_windows(args.genome, window=WINDOW)
    genes = gene_density(args.gff, nuc, window=WINDOW)
    busco = busco_orthologs(args.busco)
    org_feats = organelle_features(args.mito_gb) + organelle_features(args.pltd_gb)

    synteny = load_self_synteny(args.synteny)
    numts = load_organelle_links(args.numt) if Path(args.numt).exists() else []
    nupts = load_organelle_links(args.nupt) if Path(args.nupt).exists() else []

    # Drop any links referencing void 
    known = set(sectors)
    synteny = [l for l in synteny if l["q_chr"] in known and l["t_chr"] in known]
    numts   = [l for l in numts   if l["org_chr"] in known and l["nuc_chr"] in known]
    nupts   = [l for l in nupts   if l["org_chr"] in known and l["nuc_chr"] in known]
    print(f"[info] {len(synteny)} synteny, {len(numts)} NUMTs, {len(nupts)} NUPTs",
          file=sys.stderr)

    chrom_colors = chromosome_palette(list(nuc))
    circos = Circos(sectors, space=3)

    for sector in circos.sectors:
        name   = sector.name
        is_org = name in org

        # Outer 
        fc = organelle_color(name) if is_org else chrom_colors[name]
        sector.axis(fc=fc, ec="black", lw=0.4)
        sector.text(name, r=110, size=10, weight="bold")

        # Gap track
        gap_tr = sector.add_track((96, 99))
        gap_tr.axis(fc="white", ec="none")
        if not is_org:
            for c, gs, ge in gaps:
                if c == name:
                    gap_tr.rect(gs, ge, fc=GAP_COLOR, ec="none")

        # Repeats placeholder
        rep_tr = sector.add_track((90, 95))
        rep_tr.axis(fc="#f5f5f5", ec="black", lw=0.2)

        # Gene density 
        gd_tr = sector.add_track((80, 89))
        gd_tr.axis(fc="white", ec="black", lw=0.2)
        if not is_org:
            vals = [v for c, _, _, v in genes if c == name]
            if vals:
                gd_tr.heatmap(vals, cmap="Greens",
                              vmin=0, vmax=max(max(vals), 1))

        # GC content
        gc_tr = sector.add_track((70, 79))
        gc_tr.axis(fc="white", ec="black", lw=0.2)
        if not is_org:
            cgc = [(s, e, v) for c, s, e, v in gc if c == name]
            if cgc:
                x = np.array([(s + e) / 2 for s, e, _ in cgc])
                y = np.array([v for _, _, v in cgc])
                ymean = float(y.mean())
                gc_tr.line(x, y, vmin=y.min(), vmax=y.max(),
                           color="#1a237e", lw=0.5)
                gc_tr.fill_between(x, y, y2=np.full_like(y, ymean),
                                   vmin=y.min(), vmax=y.max(),
                                   fc="#5c6bc0", alpha=0.4)

        # BUSCO
        bu_tr = sector.add_track((60, 69))
        bu_tr.axis(fc="white", ec="black", lw=0.2)
        if not is_org:
            cbu = [(pos, st) for c, pos, st in busco if c == name]
            for status, color in BUSCO_COLORS.items():
                xs = [pos for pos, s in cbu if s == status]
                if xs:
                    ys = [BUSCO_Y[status]] * len(xs)
                    bu_tr.scatter(xs, ys, vmin=0, vmax=1,
                                  color=color, s=4, marker=".", alpha=0.8)

        # Organelle gene track
        if is_org:
            of_tr = sector.add_track((78, 90))
            of_tr.axis(fc="white", ec="black", lw=0.2)
            for c, gname, fs, fe, strand, ftype in org_feats:
                if c != name:
                    continue
                of_tr.rect(fs * ORG_SCALE, fe * ORG_SCALE,
                           fc=ORG_FEAT_COLORS.get(ftype, "#888"), ec="none")

    # Self-synteny links
    for l in synteny:
        color = SYNTENY_FWD if l["strand"] == "+" else SYNTENY_INV
        alpha = min(0.2 + l["n_genes"] / 50.0, 0.6)
        circos.link(
            (l["q_chr"], l["q_start"], l["q_end"]),
            (l["t_chr"], l["t_start"], l["t_end"]),
            color=color, alpha=alpha, lw=0.4,
        )

    # NUMT/NUPT links
    for l in numts:
        circos.link(
            (l["org_chr"], l["org_start"] * ORG_SCALE, l["org_end"] * ORG_SCALE),
            (l["nuc_chr"], l["nuc_start"], l["nuc_end"]),
            color=NUMT_COLOR, alpha=0.6, lw=0.5,
        )
    for l in nupts:
        circos.link(
            (l["org_chr"], l["org_start"] * ORG_SCALE, l["org_end"] * ORG_SCALE),
            (l["nuc_chr"], l["nuc_start"], l["nuc_end"]),
            color=NUPT_COLOR, alpha=0.6, lw=0.5,
        )

    fig = circos.plotfig(figsize=(14, 14), dpi=150)

    legend = [
        Patch(fc=GAP_COLOR,              label="Assembly gap (Ns)"),
        Patch(fc="#bdbdbd",              label="Repeats (placeholder)"),
        Patch(fc=cm.Greens(0.7),         label="Gene density"),
        Patch(fc="#5c6bc0",              label="GC content"),
        Patch(fc=BUSCO_COLORS["Complete"],   label="BUSCO Complete"),
        Patch(fc=BUSCO_COLORS["Duplicated"], label="BUSCO Duplicated"),
        Patch(fc=BUSCO_COLORS["Fragmented"], label="BUSCO Fragmented"),
        Patch(fc=SYNTENY_FWD,            label="Self-synteny (forward)"),
        Patch(fc=SYNTENY_INV,            label="Self-synteny (inverted)"),
        Patch(fc=NUMT_COLOR,             label="NUMT (mito insertion)"),
        Patch(fc=NUPT_COLOR,             label="NUPT (plastid insertion)"),
    ]
    fig.legend(
        handles=legend, loc="center", bbox_to_anchor=(0.5, 0.5),
        fontsize=8, frameon=True,
        title=f"Lolium multiflorum var Tremona\n(organelles scaled {ORG_SCALE}x)",
        title_fontsize=10,
    )

    fig.savefig(args.output, bbox_inches="tight")
    if args.output.endswith(".pdf"):
        png = args.output[:-4] + ".png"
        fig.savefig(png, bbox_inches="tight", dpi=200)
        print(f"[info] wrote {png}", file=sys.stderr)
    print(f"[info] wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--genome",       required=True)
    p.add_argument("--fai",          required=True)
    p.add_argument("--gff",          required=True)
    p.add_argument("--busco",        required=True)
    p.add_argument("--synteny",      required=True)
    p.add_argument("--numt",         required=True)
    p.add_argument("--nupt",         required=True)
    p.add_argument("--mito-fasta",   required=True, dest="mito_fasta")
    p.add_argument("--pltd-fasta",   required=True, dest="pltd_fasta")
    p.add_argument("--mito-gb",      required=True, dest="mito_gb")
    p.add_argument("--pltd-gb",      required=True, dest="pltd_gb")
    p.add_argument("--output",       default="results/data_circo/circos.pdf")
    args = p.parse_args()
    build_plot(args)