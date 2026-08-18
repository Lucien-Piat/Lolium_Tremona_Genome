import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

THRESHOLD = 0.20


ALLELIC_MEAN = 0.03
WGD_CENTRE, WGD_SD = 0.77, 0.09
ARTEFACT_FRACTION = 0.5

bp = pd.read_csv(snakemake.input.bp, sep="\t")
blocks = pd.read_csv(snakemake.input.blocks, sep="\t", low_memory=False)

het = bp[bp.experiment == "ks"].copy()
het["condition"] = het.condition.astype(float)
kb = blocks[blocks.experiment == "ks"].copy()
kb["condition"] = kb.condition.astype(float)

with open(snakemake.input.fai) as fh:
    asm_bp = sum(int(l.split("\t")[1]) for l in fh if l.strip())

floor_bp = 0
try:
    with open(snakemake.input.bed) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) >= 3:
                floor_bp += int(p[2]) - int(p[1])
except OSError:
    pass

plt.style.use("seaborn-v0_8-white")
plt.rcParams.update({
    "axes.spines.top": False, "axes.grid": False,
    "font.size": 11, "axes.labelsize": 11,
    "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
    "legend.fontsize": 9.5, "legend.frameon": False,
    "axes.labelcolor": "#212121", "text.color": "#212121",
    "savefig.dpi": 400, "savefig.bbox": "tight",
})
CORRECT, MISSED, WRONG, PALE, MUTED = ("#9ecae1", "#8a8a8a", "#081d58",
                                       "#c4c4c4", "#616161")


RECENT, ANCIENT = "#f57c00", "#2e7d32"

conds = sorted(het.condition.unique())


het["missed_bp"] = het.should_remove_bp - het.correct_bp
het["over_bp"] = het.boundary_bp + het.wrong_bp + het.background_bp


def summarise(col):
    g = het.groupby("condition")[col]
    m, s, n = g.mean(), g.std(), g.count()
    return m / 1e6, (1.96 * s / np.sqrt(n.clip(lower=1)) / 1e6).fillna(0.0)


a_correct, a_ec = summarise("correct_bp")
a_missed, a_em = summarise("missed_bp")
a_over, a_eo = summarise("over_bp")


grid = np.sort(kb.injected_ks.unique())
bp_at = kb.groupby("injected_ks")["copy_bp"].mean().reindex(grid).to_numpy()
n_per_run = len(kb) / kb.rep.nunique() / kb.condition.nunique()

density = (ARTEFACT_FRACTION * np.exp(-grid / ALLELIC_MEAN)
           + (1 - ARTEFACT_FRACTION)
           * np.exp(-0.5 * ((grid - WGD_CENTRE) / WGD_SD) ** 2))


# Weight by density times the interval each grid point stands for: the
# ladder is 40x finer near zero than in its tail, and the density alone
# would give the Ks = 0 point 18% of the total weight.
edges = np.r_[grid[0], (grid[1:] + grid[:-1]) / 2, grid[-1]]
width = np.diff(edges)
width[0] = (grid[1] - grid[0]) / 2
width[-1] = (grid[-1] - grid[-2]) / 2
w = density * width
w = w / w.sum()
below = grid < THRESHOLD

b_correct, b_missed, b_over = [], [], []
for h in conds:


    per_rep = (kb[kb.condition == h]
               .groupby(["rep", "injected_ks"])["detected"].mean()
               .unstack().reindex(columns=grid))
    p = per_rep.to_numpy()
    scale = n_per_run * w * bp_at / 1e6
    b_correct.append((p * scale * below).sum(axis=1))
    b_missed.append(((1 - p) * scale * below).sum(axis=1))
    b_over.append((p * scale * ~below).sum(axis=1))


def stats(rows):
    m = np.array([r.mean() for r in rows])
    e = np.array([1.96 * r.std(ddof=1) / np.sqrt(len(r)) for r in rows])
    return m, np.nan_to_num(e)


b_c, b_ec = stats(b_correct)
b_m, b_em = stats(b_missed)
b_o, b_eo = stats(b_over)


fig = plt.figure(figsize=(9.6, 6.4))
gs = fig.add_gridspec(2, 2, width_ratios=[2.1, 1.0], hspace=0.30, wspace=0.30)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
ax_ad = fig.add_subplot(gs[0, 1])


ax_bd = fig.add_subplot(gs[1, 1], sharex=ax_ad)

x = np.arange(len(conds))
wd = 0.27
bar_kw = dict(edgecolor="white", lw=1.4, zorder=3,
              error_kw=dict(ecolor="#4a4a4a", lw=1.0, capsize=2.5))

for ax, (c, ec), (m, em), (o, eo) in (
        (ax_a, (a_correct, a_ec), (a_missed, a_em), (a_over, a_eo)),
        (ax_b, (b_c, b_ec), (b_m, b_em), (b_o, b_eo))):
    ax.bar(x - wd, c, width=wd, yerr=ec, color=CORRECT,
           label="correctly purged", **bar_kw)
    ax.bar(x, m, width=wd, yerr=em, color=MISSED,
           label="missed (should have gone)", **bar_kw)
    ax.bar(x + wd, o, width=wd, yerr=eo, color=WRONG,
           label="over-purged (should have stayed)", **bar_kw)
    if floor_bp > 0:
        ax.axhline(floor_bp / 1e6, color=PALE, lw=1.2, ls="--", zorder=1)
    ax.set_xlim(-0.65, len(x) - 0.35)
    ax.set_ylabel("Sequence (Mb)")

top = float(max(np.max(a_correct + a_ec), np.max(a_missed + a_em),
                np.max(a_over + a_eo),
                np.max(b_c + b_ec), np.max(b_m + b_em),
                np.max(b_o + b_eo))) * 1.08
ax_a.set_ylim(0, top)
ax_b.set_ylim(0, top)
ax_a.set_title("(a)  Near-uniform $K_s$ ladder", loc="left",
               fontweight="bold", fontsize=12, pad=6)
ax_b.set_title("(b)  Bimodal $K_s$ distribution", loc="left",
               fontweight="bold", fontsize=12, pad=6)


uniform_counts = (kb[kb.condition == conds[0]]
                  .groupby("injected_ks").size().reindex(grid).to_numpy()
                  / kb.rep.nunique())
weighted_counts = n_per_run * w


BIN = 0.025
edges_h = np.arange(0.0, 1.3 + BIN, BIN)
centres = (edges_h[:-1] + edges_h[1:]) / 2

for ax, counts in ((ax_ad, uniform_counts), (ax_bd, weighted_counts)):
    binned, _ = np.histogram(grid, bins=edges_h, weights=counts)
    recent = np.where(edges_h[:-1] < THRESHOLD, binned, 0.0)
    ancient = np.where(edges_h[:-1] >= THRESHOLD, binned, 0.0)
    ax.bar(centres, recent, width=BIN, color=RECENT, align="center",
           edgecolor="none", zorder=3)
    ax.bar(centres, ancient, width=BIN, color=ANCIENT, align="center",
           edgecolor="none", zorder=3)
    ax.axvline(THRESHOLD, color=MUTED, lw=1.0, ls="--", zorder=1)
    ax.set_xlim(0, 1.3)
    ax.set_ylim(0, binned.max() * 1.18)
    ax.set_ylabel("Blocks per run")

ax_bd.set_xlabel("Injected $K_s$")
plt.setp(ax_ad.get_xticklabels(), visible=False)


ax_b.set_xticks(x)
ax_b.set_xticklabels(["Base" if c == 0 else f"+{c * 100:g}%" for c in conds])
ax_b.set_xlabel("Heterozygosity of the assembly")
plt.setp(ax_a.get_xticklabels(), visible=False)
ax_a.legend(loc="lower left", bbox_to_anchor=(0.0, 1.16), ncol=3,
            handletextpad=0.5, borderpad=0.0, columnspacing=1.6, fontsize=9.5)

for ax in (ax_a, ax_b):
    a2 = ax.twinx()
    a2.set_ylim(0, 100.0 * top * 1e6 / asm_bp)
    a2.set_ylabel(f"% of {asm_bp / 1e6:.0f} Mb")
    a2.spines["top"].set_visible(False)

fig.savefig(snakemake.output.png)
fig.savefig(snakemake.output.pdf)
plt.close(fig)

frac_boundary = w[(grid >= 0.20) & (grid <= 0.25)].sum()
for i, c in enumerate(conds):
    lab = "Base" if c == 0 else f"+{c * 100:g}%"
