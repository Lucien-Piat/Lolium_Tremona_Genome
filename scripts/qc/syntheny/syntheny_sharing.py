"""Cross-species sharing of self-synteny (duplication) blocks, via orthogroups.

Views:
1. Tremona-referenced (panels A, B, C, UpSet): each Tremona block scored
   present/absent per species. Collapse-vs-real readout for Tremona.
2. Symmetric consensus (panels D, E, networks F, G): blocks from ALL assemblies
   pooled and clustered into consensus blocks by shared orthogroup pairs, each
   with a Ks = median of member blocks (per-assembly blocks_ks.tsv from
   synteny_ks_all.py).

Colour code: red = Tremona-only, blue = shared within Lolium, green = with outgroups.
In panel A the same colours are used, with a diverging layout: intra-chromosomal
blocks above the zero line, inter-chromosomal blocks below it.
"""
import argparse
import os
import re
import sys
from collections import Counter

import matplotlib # type: ignore
matplotlib.use("Agg") # type: ignore
import matplotlib.pyplot as plt # type: ignore
import matplotlib.transforms as mtransforms # type: ignore
from matplotlib.patches import Patch # type: ignore
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
NODE_SIZE = 1000          # fixed network node size (no longer encodes a value)

BIN_WIDTH = 0.03
HEADER_RE = re.compile(r"## Alignment\s+(\d+):.*?(\S+)&(\S+)\s+(plus|minus)")

# clean display names; binomials are italicised via font style (not mathtext)
# so that font weight (bold) still applies consistently, e.g. in the networks.
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
CAT_COLOR = {"private": "#c62828", "lolium": "#1976d2", "ancient": "#2e7d32"}
BUSCO_COLOR = "#ff7f0e"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.labelsize": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 9,
})




def _net_color(sp):
    if sp == REFERENCE:
        return CAT_COLOR["private"]
    if sp in OUTGROUPS:
        return CAT_COLOR["ancient"]
    return CAT_COLOR["lolium"]


def _groups():
    outg = [s for s in OUTGROUPS if s in DATASETS]
    lol  = [s for s in DATASETS if s not in outg and s != REFERENCE]
    return outg, lol


# ----------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------
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
    ax.text(-0.07, 1.08, letter, transform=ax.transAxes, fontsize=16,
            fontweight="bold", va="bottom", ha="right")


def _draw_collapse_hist(ax, df, xmax=1.2):
    """Diverging histogram: intra-chromosomal blocks above the zero line,
    inter-chromosomal blocks below it. Fill colour still encodes category."""
    bins    = np.arange(0, xmax + BIN_WIDTH, BIN_WIDTH)
    centres = (bins[:-1] + bins[1:]) / 2

    intra = (df["q_chr"] == df["t_chr"])

    top = np.zeros(len(centres))   # running height of the intra stack (positive)
    bot = np.zeros(len(centres))   # running height of the inter stack (drawn negative)
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
    # show absolute counts on both halves of the y axis
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: str(int(abs(v)))))
    # label which half is which
    # Intra-chr annotations (Top)
    ax.text(0.48, 0.95, "■", transform=ax.transAxes, ha="right", va="top",
            fontsize=12, color="#9c27b0")
    ax.text(0.49, 0.95, "Intra-chr", transform=ax.transAxes, ha="left", va="top",
            fontsize=12, color="#555")

    # Inter-chr annotations (Bottom)
    ax.text(0.48, 0.05, "■", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=12, color="#1976d2")
    ax.text(0.49, 0.05, "Inter-chr", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=12, color="#555")


def _draw_busco(ax, ks_tsv, busco_path, xmax=1.2, ks_split=0.2):
    try:
        bt = pd.read_csv(ks_tsv, sep="\t")
    except Exception:
        bt = pd.DataFrame()
    have_coords = {"q_chr", "q_start", "q_end", "t_chr", "t_start",
                   "t_end", "median_ks"}.issubset(bt.columns)
    ax.set_ylabel("Tremona duplicated\nBUSCOs %", size='small')
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
    seen, cum = set(), []
    for k in range(len(bins) - 1):
        lo, hi = bins[k], bins[k + 1]
        for r in recs:
            mk = r.get("median_ks", np.nan)
            if pd.isna(mk) or not (lo <= mk < hi):
                continue
            for chrom, s, e in ((r["q_chr"], r["q_start"], r["q_end"]),
                                (r["t_chr"], r["t_start"], r["t_end"])):
                if pd.isna(s) or s == -1:
                    continue
                over = dup[(dup["Sequence"] == chrom) &
                           (dup["Gene_Start"] <= e) & (dup["Gene_End"] >= s)]
                seen.update(over["Busco_id"].tolist())
        cum.append(len(seen))
    pct = (np.array(cum, dtype=float) / total * 100) if total else np.zeros(len(cum))

    ax.fill_between(centres, 0, pct, color=BUSCO_COLOR, alpha=0.2)
    ax.plot(centres, pct, color=BUSCO_COLOR, lw=2)
    ax.axvline(ks_split, color="#222", ls="--", lw=1)


def _draw_dot(ax, df, xmax=1.2):
    breadth = df[DATASETS].astype(int).sum(axis=1).values
    jitter  = (np.random.RandomState(0).rand(len(df)) - 0.5) * 0.3
    ax.scatter(df["ks"].values, breadth + jitter,
               c=[CAT_COLOR[c] for c in df["category"]],
               s=18, alpha=0.8, edgecolor="none")
    ax.set_ylabel("Species sharing")
    ax.set_xlim(0, xmax)
    ax.set_ylim(0.5, len(DATASETS) + 0.5)
    ax.set_yticks(range(1, len(DATASETS) + 1))


def _draw_consensus_heat(ax, df_con, mode, cmap_name, ylabel, xmax, ks_split):
    rows, M = _consensus_counts(df_con, mode, xmax)

    # sort rows by their Sigma total, independently per panel (largest on top)
    order = np.argsort(M.sum(axis=1), kind="stable")
    rows = [rows[i] for i in order]
    M = M[order]

    nb = M.shape[1]
    binw = xmax / nb
    vmax = max(float(M.max()), 1.0)

    cm = matplotlib.colormaps[cmap_name].copy()
    cm.set_bad("white")                       # 0 cells are pure white
    Mm = np.ma.masked_where(M == 0, M)
    ax.imshow(Mm, aspect="auto", origin="lower", cmap=cm, vmin=0, vmax=vmax,
              extent=[0, xmax, 0, len(rows)])

    # number in each non-empty cell
    for j in range(len(rows)):
        for k in range(nb):
            v = M[j, k]
            if v <= 0:
                continue
            tc = "white" if v / vmax > 0.55 else "#222"
            ax.text((k + 0.5) * binw, j + 0.5, str(int(v)),
                    ha="center", va="center", fontsize=7, color=tc)

    # row totals to the right (Sigma column)
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    totals = M.sum(axis=1).astype(int)
    for j in range(len(rows)):
        ax.text(1.015, j + 0.5, str(int(totals[j])), transform=trans,
                ha="left", va="center", fontsize=9, fontweight="bold",
                color="#222", clip_on=False)
    ax.text(1.015, len(rows) + 0.12, "Σ", transform=trans, ha="left",
            va="bottom", fontsize=10, fontweight="bold", color="#222", clip_on=False)

    ax.set_yticks(np.arange(len(rows)) + 0.5)
    ax.set_yticklabels([disp(r) for r in rows])
    for t, r in zip(ax.get_yticklabels(), rows):
        t.set_fontstyle(disp_style(r))
    ax.set_ylabel(ylabel, size = "small")
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
            ax.text(m[0], m[1], str(w), fontsize=7.5, ha="center", va="center",
                    color="#333", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
        for sp in order:
            x, y = pos[sp]
            ax.scatter(x, y, s=NODE_SIZE, color=_net_color(sp), zorder=3,
                       edgecolor="white", linewidth=1.8)
        for sp in order:
            x, y = pos[sp]
            a = ang_of[sp]
            ha = "left" if np.cos(a) > 0.15 else ("right" if np.cos(a) < -0.15 else "center")
            va = "bottom" if np.sin(a) > 0.15 else ("top" if np.sin(a) < -0.15 else "center")
            ax.text(x * 1.16, y * 1.16, disp(sp), ha=ha, va=va, fontsize=9,
                    fontweight="bold", fontstyle=disp_style(sp), color="#222")

    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(-1.6, 1.7)
    ax.set_aspect("auto")
    ax.axis("off")
    ax.set_title(title, fontsize=12, fontweight="bold")


# ----------------------------------------------------------------------
# UpSet (Tremona-referenced)
# ----------------------------------------------------------------------
def fig_upset(df_sub, outpath, title):
    if df_sub.empty:
        print(f"[warn] no blocks for {title}", file=sys.stderr)
        return
    outg, _ = _groups()
    sigs = Counter()
    for _, row in df_sub.iterrows():
        sig = tuple(s for s in DATASETS if row[s])
        if sig:
            sigs[sig] += 1
    items = sorted(sigs.items(), key=lambda kv: -kv[1])
    xs = np.arange(len(items))
    bar_colors = ["#1976d2" if any(s in outg for s in sig) else "#9e9e9e"
                  for sig, _ in items]

    fig = plt.figure(figsize=(max(6, len(items) * 0.6), 5.5))
    gs  = fig.add_gridspec(2, 1, height_ratios=[3, 2.2], hspace=0.05)
    axb = fig.add_subplot(gs[0])
    axm = fig.add_subplot(gs[1], sharex=axb)
    axb.bar(xs, [c for _, c in items], color=bar_colors, width=0.7)
    for x, (_, c) in zip(xs, items):
        axb.text(x, c, str(c), ha="center", va="bottom", fontsize=9)
    axb.set_ylabel("blocks")
    axb.set_title(title)
    axb.set_xticks([])
    axb.legend(handles=[Patch(color="#1976d2", label="reaches an outgroup"),
                        Patch(color="#9e9e9e", label="Lolium only")],
               loc="upper right", fontsize=9)
    for yi, sp in enumerate(DATASETS):
        y = len(DATASETS) - 1 - yi
        for x, (sig, _) in zip(xs, items):
            axm.scatter(x, y, s=110, color=("#222" if sp in sig else "#e0e0e0"), zorder=3)
    for x, (sig, _) in zip(xs, items):
        ys = [len(DATASETS) - 1 - DATASETS.index(s) for s in sig]
        if len(ys) > 1:
            axm.plot([x, x], [min(ys), max(ys)], color="#222", lw=2, zorder=2)
    axm.set_yticks(range(len(DATASETS)))
    axm.set_yticklabels([disp(s) for s in DATASETS[::-1]], fontsize=10)
    for t, s in zip(axm.get_yticklabels(), DATASETS[::-1]):
        t.set_fontstyle(disp_style(s))
    axm.set_xticks([])
    axm.set_ylim(-0.5, len(DATASETS) - 0.5)
    axm.spines["left"].set_visible(False)
    axm.spines["bottom"].set_visible(False)
    fig.savefig(outpath, bbox_inches="tight")
    fig.savefig(outpath.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------
# composite
# ----------------------------------------------------------------------
def fig_overview(df_ref, df_con, ks_split, outpath, ks_tsv, busco_path, xmax=1.2):
    df_ref = classify_blocks(df_ref)
    recent = df_con[df_con["ks"] <  ks_split]
    old    = df_con[df_con["ks"] >= ks_split]

    fig = plt.figure(figsize=(16, 10))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.14)
    left  = outer[0].subgridspec(5, 1, height_ratios=[1.4, 0.7, 0.9, 0.55, 0.55],
                                 hspace=0.22)
    right = outer[1].subgridspec(2, 1, hspace=0.04)

    axA  = fig.add_subplot(left[0])
    axB  = fig.add_subplot(left[1], sharex=axA)   # BUSCO
    axC  = fig.add_subplot(left[2], sharex=axA)   # dots
    axDd = fig.add_subplot(left[3], sharex=axA)   # duplicated heatmap
    axDm = fig.add_subplot(left[4], sharex=axA)   # missing heatmap

    _draw_collapse_hist(axA, df_ref, xmax)
    _draw_busco(axB, ks_tsv, busco_path, xmax, ks_split)
    _draw_dot(axC, df_ref, xmax)
    _draw_consensus_heat(axDd, df_con, "dup", "Reds", "Lineage-specific", xmax, ks_split)
    _draw_consensus_heat(axDm, df_con, "missing", "Greens", "Lost", xmax, ks_split)

    for ax in (axA, axC):
        ax.axvline(ks_split, color="#222", ls="--", lw=1)
    axA.text(ks_split, axA.get_ylim()[1] * 0.98, " old", fontsize=8,
             va="top", ha="left", color="#555")
    axA.text(ks_split, axA.get_ylim()[1] * 0.98, "recent ", fontsize=8,
             va="top", ha="right", color="#555")
    axA.legend(handles=[Patch(color=CAT_COLOR[c], label=CAT_LABEL[c])
                        for c in CAT_ORDER], loc="upper right", fontsize=9)

    for ax in (axA, axB, axC, axDd):
        ax.tick_params(labelbottom=False)
    axDm.set_xlabel("Median Ks (NG86)")

    axF = fig.add_subplot(right[0])
    axG = fig.add_subplot(right[1])
    _draw_network(axF, recent, f"Recent  (Ks < {ks_split})")
    _draw_network(axG, old,    f"Old  (Ks >= {ks_split})")

    for ax, L in [(axA, "A.        "), (axB, "B.        "), (axC, "C.       "),
                  (axDd, "D.        "), (axDm, "E.        "), (axF, "F."), (axG, "G.")]:
        _panel_letter(ax, L)

    fig.savefig(outpath, bbox_inches="tight")
    fig.savefig(outpath.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_cytoscape(df, outdir):
    with open(os.path.join(outdir, "cytoscape_edges.tsv"), "w") as e:
        e.write("block\tspecies\n")
        for bid, row in df.iterrows():
            for sp in DATASETS:
                if row[sp]:
                    e.write(f"B{bid}\t{sp}\n")
    with open(os.path.join(outdir, "cytoscape_nodes.tsv"), "w") as n:
        n.write("node\ttype\tks\tbreadth\n")
        for sp in DATASETS:
            n.write(f"{sp}\tspecies\tNA\tNA\n")
        for bid, row in df.iterrows():
            breadth = int(sum(bool(row[sp]) for sp in DATASETS))
            ks = row["ks"] if not pd.isna(row["ks"]) else float("nan")
            n.write(f"B{bid}\tblock\t{ks:.4f}\t{breadth}\n")


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

    recent = df_ref[df_ref["ks"] <  args.ks_split]
    old    = df_ref[df_ref["ks"] >= args.ks_split]
    fig_upset(recent, os.path.join(args.outdir, "upset_recent.pdf"),
              f"Recent blocks (Ks < {args.ks_split})")
    fig_upset(old, os.path.join(args.outdir, "upset_old.pdf"),
              f"Old blocks (Ks >= {args.ks_split})")
    fig_overview(df_ref, df_con, args.ks_split,
                 os.path.join(args.outdir, "dupshare_overview.pdf"),
                 args.ks_tsv, args.busco)
    write_cytoscape(df_con, args.outdir)
    print(f"[info] done -> {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Orthogroup-based sharing of synteny blocks.")
    p.add_argument("--orthogroups", required=True, help="OrthoFinder Orthogroups.tsv")
    p.add_argument("--syn-base", default="results/synteny")
    p.add_argument("--ks-tsv", default="results/synteny/tremona/synteny_ks.tsv")
    p.add_argument("--busco",
                   default="reference_data/lmultiflorum.tremona_full_table_busco_format.tsv")
    p.add_argument("--outdir", default="results/dupshare")
    p.add_argument("--reference", default="tremona")
    p.add_argument("--ks-split", type=float, default=0.2)
    p.add_argument("--min-shared", type=int, default=2)
    p.add_argument("--match-frac", type=float, default=0.3)
    main(p.parse_args())