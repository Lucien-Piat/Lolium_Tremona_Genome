"""
Cross-species sharing of self-synteny (duplication) blocks, via orthogroups.
"""
import argparse
import os
import re
import sys

import matplotlib # type: ignore
matplotlib.use("Agg") # type: ignore
import matplotlib.pyplot as plt # type: ignore
import matplotlib.transforms as mtransforms # type: ignore
from matplotlib.lines import Line2D # type: ignore
from matplotlib.patches import Rectangle # type: ignore
from matplotlib.ticker import FuncFormatter # type: ignore
import numpy as np # type: ignore
import pandas as pd # type: ignore

DATASETS  = ["tremona", "rabiosa", "sikem", "perenne", "paraquat",
             "brachypodium", "oryza"]
REFERENCE = "tremona"
OUTGROUPS = ["brachypodium", "oryza"]

NET_ORDER = ["brachypodium", "oryza", "perenne", "tremona", "paraquat", 
             "sikem", "rabiosa"]

HEAT_ROWS = ["perenne", "paraquat", "sikem", "rabiosa", "tremona"]
HEAT_BIN  = 0.03
NODE_SIZE = 650  

BIN_WIDTH = 0.03
HEADER_RE = re.compile(r"## Alignment\s+(\d+):.*?(\S+)&(\S+)\s+(plus|minus)")

DISPLAY = {
    "tremona": "Tremona", "rabiosa": "Rabiosa", "sikem": "Sikem",
    "paraquat": "Brunharo", "perenne": "L. perenne",
    "brachypodium": "B. distachyon", "oryza": "O. sativa ",
}

ITALIC = {"perenne", "brachypodium", "oryza"}
def disp(sp):
    return DISPLAY.get(sp, sp)
def disp_style(sp):
    return "italic" if sp in ITALIC else "normal"

CAT_ORDER = ["private", "lolium", "ancient"]
CAT_LABEL = {"private": "Tremona-only",
             "lolium":  "Shared within Lolium only",
             "ancient": "Shared with outgroups"}


CAT_MARKER = {"private": "o", "lolium": "s", "ancient": "^"}
CAT_COLOR = {"private": "#e4000bff", "lolium": "#1976d2", "ancient": "#2e7d32"}
INTRA_COLOR = "#000000"
INTER_COLOR = "#000000"
BUSCO_COLOR = "#000000"

# Increased all base font sizes for readability from afar
plt.rcParams.update({
    "font.family": "DejaVu Sans", 
    "font.size": 16, 
    "axes.titlesize": 18,
    "axes.labelsize": 16, 
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.spines.top": False, 
    "axes.spines.right": False,
    "legend.frameon": False, 
    "legend.fontsize": 14,
    "hatch.linewidth": 0.5,
})

def _net_color(sp):
    if sp == REFERENCE:
        return CAT_COLOR["private"]
    if sp in OUTGROUPS:
        return CAT_COLOR["ancient"]
    return CAT_COLOR["lolium"]


def _net_marker(sp):
    if sp == REFERENCE:
        return CAT_MARKER["private"]
    if sp in OUTGROUPS:
        return CAT_MARKER["ancient"]
    return CAT_MARKER["lolium"]


def _groups():
    outg = [s for s in OUTGROUPS if s in DATASETS]
    lol  = [s for s in DATASETS if s not in outg and s != REFERENCE]
    return outg, lol


def load_orthogroups(path):
    og_of = {}
    with open(path) as fh:
        species = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            og = cells[0]
            for sp, blob in zip(species, cells[1:]):
                if not blob:
                    continue
                for gene in blob.split(", "):
                    gene = gene.strip()
                    if gene:
                        og_of[(sp, gene)] = og
    return og_of


def parse_collinearity(path):
    blocks, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("## Alignment"):
                m = HEADER_RE.search(line)
                if m:
                    cur = int(m.group(1))
                    blocks[cur] = {"pairs": [], "q_chr": m.group(2), "t_chr": m.group(3)}
            elif cur is not None and line and not line.startswith("#"):
                p = line.split("\t")
                if len(p) >= 3:
                    blocks[cur]["pairs"].append((p[1].strip(), p[2].strip()))
    return blocks


def block_ogpairs(blocks, species, og_of):
    out = {}
    for bid, b in blocks.items():
        pairs = set()
        for g1, g2 in b["pairs"]:
            o1 = og_of.get((species, g1))
            o2 = og_of.get((species, g2))
            if o1 and o2:
                pairs.add(tuple(sorted((o1, o2))))
        if pairs:
            out[bid] = {"ogpairs": pairs, "q_chr": b["q_chr"], "t_chr": b["t_chr"]}
    return out


def load_block_ks(syn_base):
    ks = {}
    for sp in DATASETS:
        d = {}
        for fn in ("blocks_ks.tsv", "synteny_ks.tsv"):
            path = os.path.join(syn_base, sp, fn)
            if os.path.isfile(path):
                try:
                    t = pd.read_csv(path, sep="\t")
                    if "block_id" in t.columns and "median_ks" in t.columns:
                        d = t.set_index("block_id")["median_ks"].to_dict()
                        break
                except Exception:
                    pass
        ks[sp] = d
        if not d:
            print(f"[warn] no per-block Ks for {sp} (run synteny_ks_all.py)", file=sys.stderr)
    return ks


# ----------------------------------------------------------------------
# view 1: Tremona-referenced presence
# ----------------------------------------------------------------------
def present_in(ref_pairs, sp_blocks, min_shared, frac):
    best = 0
    for b in sp_blocks.values():
        inter = len(ref_pairs & b["ogpairs"])
        if inter > best:
            best = inter
    return best >= min_shared and best / len(ref_pairs) >= frac


def build_presence(syn_base, ks_tsv, og_of, reference, min_shared, frac):
    sp_blocks = {}
    for sp in DATASETS:
        col = os.path.join(syn_base, sp, f"{sp}.collinearity")
        sp_blocks[sp] = block_ogpairs(parse_collinearity(col), sp, og_of)
        print(f"[info] {sp}: {len(sp_blocks[sp])} blocks with OG-pairs", file=sys.stderr)

    ks = pd.read_csv(ks_tsv, sep="\t").set_index("block_id")["median_ks"].to_dict()
    rows = []
    for bid, b in sp_blocks[reference].items():
        if bid not in ks:
            continue
        rec = {"block": bid, "ks": ks[bid], "q_chr": b["q_chr"], "t_chr": b["t_chr"]}
        for sp in DATASETS:
            rec[sp] = present_in(b["ogpairs"], sp_blocks[sp], min_shared, frac)
        rows.append(rec)
    df = pd.DataFrame(rows).set_index("block").sort_values("ks")
    print(f"[info] {len(df)} Tremona reference blocks scored", file=sys.stderr)
    return df, sp_blocks


def classify_blocks(df):
    outg, lol = _groups()
    def cat(row):
        if any(row[s] for s in outg):
            return "ancient"
        if any(row[s] for s in lol):
            return "lolium"
        return "private"
    out = df.copy()
    out["category"] = out.apply(cat, axis=1)
    return out


# ----------------------------------------------------------------------
# view 2: symmetric consensus blocks
# ----------------------------------------------------------------------
def build_consensus(sp_blocks, syn_base, min_shared, frac):
    blk_ks = load_block_ks(syn_base)

    items = []
    for sp in DATASETS:
        for bid, b in sp_blocks[sp].items():
            items.append((sp, b["ogpairs"], bid))
    n = len(items)

    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        si = items[i][1]
        li = len(si)
        for j in range(i + 1, n):
            sj = items[j][1]
            inter = len(si & sj)
            if inter >= min_shared and inter / min(li, len(sj)) >= frac:
                union(i, j)

    clusters = {}
    for idx in range(n):
        clusters.setdefault(find(idx), []).append(idx)

    rows, n_no_ks = [], 0
    for cid, members in clusters.items():
        present = {sp: False for sp in DATASETS}
        kvals = []
        for m in members:
            sp, _, bid = items[m]
            present[sp] = True
            v = blk_ks.get(sp, {}).get(bid)
            if v is not None and not pd.isna(v):
                kvals.append(float(v))
        ks = float(np.median(kvals)) if kvals else np.nan
        if not kvals:
            n_no_ks += 1
        rec = {"cluster": cid, "ks": ks}
        rec.update(present)
        rows.append(rec)

    df = pd.DataFrame(rows).set_index("cluster").sort_values("ks")
    print(f"[info] {len(df)} consensus blocks ({n_no_ks} without any Ks)", file=sys.stderr)
    return df


def _consensus_counts(df, mode, xmax):
    rows = [s for s in HEAT_ROWS if s in DATASETS]
    outg, lol = _groups()
    both_out = np.ones(len(df), dtype=bool)
    any_out  = np.zeros(len(df), dtype=bool)
    for s in outg:
        both_out &= df[s].values.astype(bool)
        any_out  |= df[s].values.astype(bool)

    bins = np.arange(0, xmax + HEAT_BIN, HEAT_BIN)
    nb   = len(bins) - 1
    ks   = df["ks"].values
    bidx = np.digitize(ks, bins) - 1

    M = np.zeros((len(rows), nb))
    for j, sp in enumerate(rows):
        pres = df[sp].values.astype(bool)
        for i in range(len(df)):
            if np.isnan(ks[i]) or not (0 <= bidx[i] < nb):
                continue
            if mode == "dup":
                if pres[i] and not any_out[i]:
                    M[j, bidx[i]] += 1
            else:
                other_lol = any(df[s].values[i] for s in lol if s != sp)
                if (not pres[i]) and both_out[i] and other_lol:
                    M[j, bidx[i]] += 1
    return rows, M


# ----------------------------------------------------------------------
# panel drawers
# ----------------------------------------------------------------------
def _panel_letter(ax, letter):
    # Increased panel letter size
    ax.text(-0.07, 1.08, letter, transform=ax.transAxes, fontsize=24,
            fontweight="bold", va="bottom", ha="right")


def _draw_collapse_hist(ax, df, xmax=1.2):
    """Diverging histogram: intra-chromosomal blocks above the zero line,
    inter-chromosomal blocks below it. Fill colour still encodes category."""
    bins    = np.arange(0, xmax + BIN_WIDTH, BIN_WIDTH)
    centres = (bins[:-1] + bins[1:]) / 2

    intra = (df["q_chr"] == df["t_chr"])

    top = np.zeros(len(centres))
    bot = np.zeros(len(centres))
    for c in CAT_ORDER:
        col = CAT_COLOR[c]
        in_counts, _  = np.histogram(df.loc[(df["category"] == c) &  intra, "ks"].values, bins=bins)
        out_counts, _ = np.histogram(df.loc[(df["category"] == c) & ~intra, "ks"].values, bins=bins)
        ax.bar(centres, in_counts, width=BIN_WIDTH * 0.95, bottom=top,
               color=col, edgecolor="white", linewidth=0.2)
        ax.bar(centres, -out_counts, width=BIN_WIDTH * 0.95, bottom=-bot,
               color=col, edgecolor="white", linewidth=0.2)
        top += in_counts
        bot += out_counts

    ax.axhline(0, color="#222", lw=0.8)
    ax.set_ylabel("Number\nof blocks")
    ax.set_xlim(0, xmax)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: str(int(abs(v)))))
    
    sw_w, sw_h = 0.035, 0.13
    ax.add_patch(Rectangle((0.48 - sw_w, 0.93 - sw_h), sw_w, sw_h,
                           transform=ax.transAxes, facecolor="none",
                           edgecolor=INTRA_COLOR, hatch="//", linewidth=0.0,
                           clip_on=False))
    # Increased font size for intra/inter legend text
    ax.text(0.3, 0.90, "Intra-chr", transform=ax.transAxes, ha="left", va="top",
            fontsize=16, color="#555")
    ax.add_patch(Rectangle((0.48 - sw_w, 0.08), sw_w, sw_h,
                           transform=ax.transAxes, facecolor="none",
                           edgecolor=INTER_COLOR, hatch="...", linewidth=0.0,
                           clip_on=False))
    ax.text(0.3, 0.1, "Inter-chr", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=16, color="#555")


def _draw_busco(ax, ks_tsv, busco_path, xmax=1.2, ks_split=0.2):
    """
    Cumulative fraction of Tremona duplicated BUSCOs covered by synteny blocks
    """
    try:
        bt = pd.read_csv(ks_tsv, sep="\t")
    except Exception:
        bt = pd.DataFrame()
    ax.set_ylabel("Tremona duplicated\nBUSCOs %\n", size='medium')
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 100)

    busco = pd.read_csv(busco_path, sep="\t", comment="#",
                        names=["Busco_id", "Status", "Sequence", "Gene_Start",
                               "Gene_End", "Strand", "Score", "Length"])
    busco["Gene_Start"] = pd.to_numeric(busco["Gene_Start"], errors="coerce")
    busco["Gene_End"]   = pd.to_numeric(busco["Gene_End"],   errors="coerce")
    dup = busco[busco["Status"] == "Duplicated"].dropna(subset=["Gene_Start", "Gene_End"])
    total = dup["Busco_id"].nunique()

    bins = np.arange(0, xmax + BIN_WIDTH, BIN_WIDTH)
    centres = (bins[:-1] + bins[1:]) / 2
    recs = bt.to_dict("records")

    seen = set()
    s_intra, s_inter = set(), set()
    cum_intra, cum_inter = [], []
    for k in range(len(bins) - 1):
        lo, hi = bins[k], bins[k + 1]
        hit_intra, hit_inter = set(), set()
        for r in recs:
            mk = r.get("median_ks", np.nan)
            if pd.isna(mk) or not (lo <= mk < hi):
                continue
            is_intra = (r.get("q_chr") == r.get("t_chr"))
            ids = set()
            for chrom, s, e in ((r["q_chr"], r["q_start"], r["q_end"]),
                                (r["t_chr"], r["t_start"], r["t_end"])):
                if pd.isna(s) or s == -1:
                    continue
                over = dup[(dup["Sequence"] == chrom) &
                           (dup["Gene_Start"] <= e) & (dup["Gene_End"] >= s)]
                ids.update(over["Busco_id"].tolist())
            if is_intra:
                hit_intra |= ids
            else:
                hit_inter |= ids

        new_intra = hit_intra - seen
        s_intra |= new_intra
        seen |= new_intra
        new_inter = hit_inter - seen
        s_inter |= new_inter
        seen |= new_inter

        cum_intra.append(len(s_intra))
        cum_inter.append(len(s_inter))

    cum_intra = np.array(cum_intra, dtype=float)
    cum_inter = np.array(cum_inter, dtype=float)
    if total:
        pct_intra = cum_intra / total * 100.0
        pct_inter = cum_inter / total * 100.0
    else:
        pct_intra = np.zeros(len(cum_intra))
        pct_inter = np.zeros(len(cum_inter))

    recent_idx = np.where(centres < ks_split)[0]
    i_split = recent_idx[-1] if len(recent_idx) else -1
    rec_intra = int(cum_intra[i_split]) if i_split >= 0 else 0
    rec_inter = int(cum_inter[i_split]) if i_split >= 0 else 0
    tot_intra = int(cum_intra[-1]) if len(cum_intra) else 0
    tot_inter = int(cum_inter[-1]) if len(cum_inter) else 0
    print(f"[info] duplicated BUSCOs covered by blocks "
          f"(of {total} duplicated in the table):", file=sys.stderr)
    print(f"[info]   recent (Ks < {ks_split}): {rec_intra + rec_inter} "
          f"(intra {rec_intra}, inter {rec_inter})", file=sys.stderr)
    print(f"[info]   old    (Ks >= {ks_split}): "
          f"{(tot_intra - rec_intra) + (tot_inter - rec_inter)} "
          f"(intra {tot_intra - rec_intra}, inter {tot_inter - rec_inter})", file=sys.stderr)
    print(f"[info]   total covered: {tot_intra + tot_inter} "
          f"(intra {tot_intra}, inter {tot_inter})", file=sys.stderr)

    ax.fill_between(centres, 0, pct_intra, facecolor="none",
                    hatch="//", edgecolor=INTRA_COLOR, linewidth=0.0)
    ax.fill_between(centres, pct_intra, pct_intra + pct_inter,
                    facecolor="none", hatch="..", edgecolor=INTER_COLOR,
                    linewidth=0.0)
    ax.plot(centres, pct_intra, color="#00000065", lw=0.4) 
    ax.plot(centres, pct_intra + pct_inter, color=BUSCO_COLOR, lw=2)  
    ax.axvline(ks_split, color="#222", ls="--", lw=1)


def _draw_consensus_heat(ax, df_con, mode, cmap_name, ylabel, xmax, ks_split):
    rows, M = _consensus_counts(df_con, mode, xmax)

    order = np.argsort(M.sum(axis=1), kind="stable")
    rows = [rows[i] for i in order]
    M = M[order]

    nb = M.shape[1]
    binw = xmax / nb
    vmax = max(float(M.max()), 1.0)

    cm = matplotlib.colormaps[cmap_name].copy()
    cm.set_bad("white")
    Mm = np.ma.masked_where(M == 0, M)
    ax.imshow(Mm, aspect="auto", origin="lower", cmap=cm, vmin=0, vmax=vmax,
              extent=[0, xmax, 0, len(rows)])

    for j in range(len(rows)):
        for k in range(nb):
            v = M[j, k]
            if v <= 0:
                continue
            tc = "white" if v / vmax > 0.55 else "#222"
            # Increased font size for heatmap cell counts
            ax.text((k + 0.5) * binw, j + 0.5, str(int(v)),
                    ha="center", va="center", fontsize=12, color=tc)

    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    totals = M.sum(axis=1).astype(int)
    for j in range(len(rows)):
        # Increased font size for row totals
        ax.text(1.015, j + 0.5, str(int(totals[j])), transform=trans,
                ha="left", va="center", fontsize=14, fontweight="bold",
                color="#222", clip_on=False)
    # Increased font size for sigma
    ax.text(1.015, len(rows) + 0.12, "Σ", transform=trans, ha="left",
            va="bottom", fontsize=16, fontweight="bold", color="#222", clip_on=False)

    ax.set_yticks(np.arange(len(rows)) + 0.5)
    ax.set_yticklabels([disp(r) for r in rows])
    for t, r in zip(ax.get_yticklabels(), rows):
        t.set_fontstyle(disp_style(r))
    ax.set_ylabel(ylabel, size="medium")
    ax.set_xlim(0, xmax)
    ax.axvline(ks_split, color="#222", ls="--", lw=1)


def _draw_network(ax, df_sub, title):
    order = NET_ORDER
    n = len(order)
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2
    pos = {sp: np.array([np.cos(a), np.sin(a)]) for sp, a in zip(order, ang)}
    ang_of = {sp: a for sp, a in zip(order, ang)}

    if len(df_sub):
        P = df_sub[order].astype(int)
        pairs = [(order[i], order[j]) for i in range(n) for j in range(i + 1, n)]
        wt = {(a, b): int((P[a] & P[b]).sum()) for a, b in pairs}
        wmax = max([1] + list(wt.values()))
        for (a, b), w in wt.items():
            if w == 0:
                continue
            p1, p2 = pos[a], pos[b]
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color="#d2d2d2",
                    lw=0.6 + 7 * w / wmax, alpha=0.85, zorder=1, solid_capstyle="round")
        for (a, b), w in wt.items():
            if w == 0:
                continue
            m = (pos[a] + pos[b]) / 2
            # Increased font size for network edge weights
            ax.text(m[0], m[1], str(w), fontsize=12, ha="center", va="center",
                    color="#333", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
        for sp in order:
            x, y = pos[sp]
            ax.scatter(x, y, s=NODE_SIZE, marker=_net_marker(sp),
                       color=_net_color(sp), zorder=3,
                       edgecolor="white", linewidth=1.8)
        for sp in order:
            x, y = pos[sp]
            a = ang_of[sp]
            ha = "left" if np.cos(a) > 0.15 else ("right" if np.cos(a) < -0.15 else "center")
            va = "bottom" if np.sin(a) > 0.15 else ("top" if np.sin(a) < -0.15 else "center")
            # Increased font size for network species labels
            ax.text(x * 1.16, y * 1.16, disp(sp), ha=ha, va=va, fontsize=14,
                    fontweight="bold", fontstyle=disp_style(sp), color="#222")

    lim = 1.65
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="datalim")  
    ax.axis("off")
    # Increased network title size
    ax.set_title(title, fontsize=18, fontweight="bold")


# ----------------------------------------------------------------------
# composite
# ----------------------------------------------------------------------
def fig_overview(df_ref, df_con, ks_split, outpath, ks_tsv, busco_path, xmax=1.2):
    df_ref = classify_blocks(df_ref)
    recent = df_con[df_con["ks"] <  ks_split]
    old    = df_con[df_con["ks"] >= ks_split]

    fig = plt.figure(figsize=(17, 10))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.14)
    
    # Removed the plot C ratio (was 0.9) to make the left column 4 plots
    left  = outer[0].subgridspec(4, 1, height_ratios=[1.4, 0.7, 0.55, 0.55],
                                 hspace=0.22)
    right = outer[1].subgridspec(2, 1, hspace=0.4)

    axA  = fig.add_subplot(left[0])
    axB  = fig.add_subplot(left[1], sharex=axA)   # BUSCO
    axDd = fig.add_subplot(left[2], sharex=axA)   # duplicated heatmap
    axDm = fig.add_subplot(left[3], sharex=axA)   # missing heatmap

    _draw_collapse_hist(axA, df_ref, xmax)
    _draw_busco(axB, ks_tsv, busco_path, xmax, ks_split)
    
    # Removed the call to `_draw_dot` for the C panel
    
    _draw_consensus_heat(axDd, df_con, "dup", "Reds", "", xmax, ks_split)
    _draw_consensus_heat(axDm, df_con, "missing", "Greens", "", xmax, ks_split)

    # Added axline directly to axA (since axC is gone)
    axA.axvline(ks_split, color="#222", ls="--", lw=1)
    
    # Increased text sizes for threshold labels
    axA.text(ks_split, axA.get_ylim()[1] * 0.98, " old", fontsize=12,
             va="top", ha="left", color="#555")
    axA.text(ks_split, axA.get_ylim()[1] * 0.98, "recent ", fontsize=12,
             va="top", ha="right", color="#555")
             
    # Increased legend marker size and label size handled by rcParams
    axA.legend(handles=[Line2D([0], [0], marker=CAT_MARKER[c], linestyle="none",
                               markerfacecolor=CAT_COLOR[c], markeredgecolor="none",
                               markersize=12, label=CAT_LABEL[c])
                        for c in CAT_ORDER], loc="upper right", fontsize=14)

    # Tick params loop minus axC
    for ax in (axA, axB, axDd):
        ax.tick_params(labelbottom=False)
    axDm.set_xlabel("Median Ks (NG86)")

    axF = fig.add_subplot(right[0])
    axG = fig.add_subplot(right[1])
    _draw_network(axF, recent, f"Recent  (Ks < {ks_split})\n")
    _draw_network(axG, old,    f"Old  (Ks >= {ks_split})\n")

    # Shifted lettering due to plot C removal
    for ax, L in [(axA, "A.        "), (axB, "B.        "),
                  (axDd, "C.        "), (axDm, "D.        "), (axF, "  E."), (axG, "  F.")]:
        _panel_letter(ax, L)

    fig.savefig(outpath, bbox_inches="tight")
    fig.savefig(outpath.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main(args):
    os.makedirs(args.outdir, exist_ok=True)
    og_of = load_orthogroups(args.orthogroups)
    print(f"[info] {len(og_of)} (species, gene) -> OG assignments", file=sys.stderr)

    df_ref, sp_blocks = build_presence(args.syn_base, args.ks_tsv, og_of,
                                       args.reference, args.min_shared, args.match_frac)
    df_ref = classify_blocks(df_ref)
    df_ref.to_csv(os.path.join(args.outdir, "presence_matrix.tsv"), sep="\t")

    df_con = build_consensus(sp_blocks, args.syn_base, args.min_shared, args.match_frac)
    df_con.to_csv(os.path.join(args.outdir, "consensus_matrix.tsv"), sep="\t")

    n_priv = int((df_ref["category"] == "private").sum())
    n_lol  = int((df_ref["category"] == "lolium").sum())
    n_anc  = int((df_ref["category"] == "ancient").sum())
    print(f"[info] Tremona categories: {n_priv} Tremona-only, {n_lol} Lolium-only, "
          f"{n_anc} reaching an outgroup", file=sys.stderr)

    fig_overview(df_ref, df_con, args.ks_split,
                 os.path.join(args.outdir, "dupshare_overview.pdf"),
                 args.ks_tsv, args.busco)
    print(f"[info] done -> {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Orthogroup-based sharing of synteny blocks.")
    p.add_argument("--orthogroups", required=True, help="OrthoFinder Orthogroups.tsv")
    p.add_argument("--syn-base", default="results/synteny")
    p.add_argument("--ks-tsv", default="results/synteny/tremona/synteny_ks.tsv")
    p.add_argument("--busco",
                   default="reference_data/lmultiflorum.tremona_full_table_busco_format.old.tsv")
    p.add_argument("--outdir", default="results/dupshare2")
    p.add_argument("--reference", default="tremona")
    p.add_argument("--ks-split", type=float, default=0.2)
    p.add_argument("--min-shared", type=int, default=2)
    p.add_argument("--match-frac", type=float, default=0.3)
    main(p.parse_args())