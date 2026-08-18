import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

THRESHOLD = float(snakemake.params.threshold)
blocks = pd.read_csv(snakemake.input.blocks, sep="\t")
scaling = pd.read_csv(snakemake.input.scaling, sep="\t")

plt.style.use("seaborn-v0_8-white")
plt.rcParams.update({
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "font.size": 11, "axes.labelsize": 11,
    "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
    "legend.fontsize": 9.5, "legend.frameon": False, "legend.handlelength": 1.5,
    "lines.linewidth": 1.4,
    "axes.labelcolor": "#212121", "text.color": "#212121",
    "savefig.dpi": 400, "savefig.bbox": "tight",
})

MUTED, PALE = "#9e9e9e", "#c4c4c4"
HETS = sorted(blocks[blocks.experiment == "ks"].condition.unique())


RAMP = plt.cm.YlGnBu(np.linspace(1.0, 0.40, len(HETS)))
HET_COLOUR = dict(zip(HETS, RAMP))
CURVE = RAMP[0]


def label(h):
    return "Base" if float(h) == 0 else f"+{float(h) * 100:g}%"


def key(ax, letter):
    ax.set_title(f"({letter})", loc="left", fontweight="bold", fontsize=13,
                 pad=6)


def curve(ax, x, y, n=400, **kw):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~np.isnan(y)
    x, y = x[ok], y[ok]
    o = np.argsort(x)
    x, y = x[o], y[o]
    if len(x) < 3:
        return ax.plot(x, y, **kw)
    xs = np.linspace(x.min(), x.max(), n)
    return ax.plot(xs, PchipInterpolator(x, y)(xs), **kw)


# Intervals are across replicate runs. Blocks within one injection share a
# genome and one self-blast, so they are not independent.
def across_reps(df, group, value):
    per_run = df.groupby(["rep", group])[value].mean().reset_index()
    g = per_run.groupby(group)[value]
    out = g.agg(["mean", "std", "count"]).reset_index()
    se = out["std"] / np.sqrt(out["count"].clip(lower=1))
    out["lo"] = (out["mean"] - 1.96 * se).clip(0, 1)
    out["hi"] = (out["mean"] + 1.96 * se).clip(0, 1)
    return out


def moving_average(y, window=3):
    y = np.asarray(y, float)
    if len(y) < window:
        return y
    pad = window // 2
    return np.convolve(np.r_[[y[0]] * pad, y, [y[-1]] * pad],
                       np.ones(window) / window, mode="valid")


fig = plt.figure(figsize=(7.09, 7.4))
gs = fig.add_gridspec(3, 100, height_ratios=[1.0, 1.0, 0.72],
                      hspace=0.62, wspace=0.0)
ax_a = fig.add_subplot(gs[0, 0:33])
ax_b = fig.add_subplot(gs[0, 44:77])
ax_key = fig.add_subplot(gs[0, 79:100]); ax_key.axis("off")
ax_c = fig.add_subplot(gs[1, 0:42])
ax_d = fig.add_subplot(gs[1, 58:100])
ax_e = fig.add_subplot(gs[2, 0:42])
ax_f = fig.add_subplot(gs[2, 58:100])


ks = blocks[blocks.experiment == "ks"]
n_blocks = ks[ks.condition == HETS[0]].groupby(["rep", "injected_ks"]).size()

ax_a.plot([0, 1.55], [0, 1.55], color=PALE, lw=1.0, ls="--", zorder=1)
for h in HETS:
    sub = ks[ks.condition == h]
    got = sub.dropna(subset=["estimated_ks"])
    est = got.groupby("injected_ks")["estimated_ks"].median()


    rate = (got.groupby("injected_ks").size()
            / sub.groupby("injected_ks").size())
    est = est[rate.reindex(est.index, fill_value=0) >= 0.5]
    curve(ax_a, est.index, est.values, color=HET_COLOUR[h], label=label(h),
          zorder=3)

    call = sub[sub.injected_ks > 0].groupby("injected_ks")["detected"].mean()
    curve(ax_b, call.index, call.values, color=HET_COLOUR[h], zorder=3)

ax_a.set_xlim(0, 1.5); ax_a.set_ylim(0, 1.55)
ax_a.set_xticks([0, 0.5, 1.0, 1.5]); ax_a.set_yticks([0, 0.5, 1.0, 1.5])
ax_a.set_ylabel("Estimated $K_s$")
key(ax_a, "a")

ax_b.axvline(THRESHOLD, color=PALE, lw=1.0, ls="--", zorder=1)
ax_b.set_xlim(0, 0.4); ax_b.set_ylim(-0.04, 1.06)
ax_b.set_xticks([0, 0.1, 0.2, 0.3, 0.4]); ax_b.set_yticks([0, 0.5, 1.0])
ax_b.set_ylabel("Blocks removed")
key(ax_b, "b")


genes = blocks[blocks.experiment == "genes"]
dep = across_reps(genes, "n_genes", "detected")
dep = dep[dep.n_genes <= 20]
ax_c.fill_between(dep.n_genes, dep.lo, dep.hi, color=CURVE, alpha=0.16, lw=0,
                  zorder=1)
curve(ax_c, dep.n_genes, dep["mean"], color=CURVE, zorder=3)
ax_c.set_xlim(0, 20.5); ax_c.set_ylim(-0.04, 1.06)
ax_c.set_xticks([0, 5, 10, 15, 20]); ax_c.set_yticks([0, 0.5, 1.0])
ax_c.set_xlabel("Genes in the block"); ax_c.set_ylabel("Recall")
key(ax_c, "c")


art = blocks[blocks.type == "artefact"]
intact = across_reps(art[art.experiment == "main"], "condition", "detected")
base_recall = float(intact["mean"].iloc[0])


LABEL_AT = {"drop": (36, 0.45), "frag": (43, 0.90)}
LABEL_OF = {"drop": "genes removed", "frag": "genes truncated"}

for mode, style in (("drop", "-"), ("frag", "--")):
    sub = art[art.experiment == f"degrade_{mode}"].copy()
    sub["condition"] = sub.condition.astype(float)
    d = across_reps(sub, "condition", "detected")
    x = np.r_[0.0, d.condition.to_numpy() * 100]
    y = np.r_[base_recall, d["mean"].to_numpy()]
    lo = np.r_[base_recall, d.lo.to_numpy()]
    hi = np.r_[base_recall, d.hi.to_numpy()]
    ax_d.fill_between(x, lo, hi, color=CURVE, alpha=0.14, lw=0, zorder=1)
    curve(ax_d, x, y, color=CURVE, ls=style, lw=1.6 if mode == "frag" else 1.4,
          zorder=3)
    ax_d.text(*LABEL_AT[mode], LABEL_OF[mode], fontsize=8.5, color=CURVE,
              ha="left", va="center")
ax_d.set_xlim(-3, 103); ax_d.set_ylim(-0.04, 1.06)
ax_d.set_xticks([0, 25, 50, 75, 100]); ax_d.set_yticks([0, 0.5, 1.0])
ax_d.set_xlabel("Annotation degraded (%)"); ax_d.set_ylabel("Recall")
key(ax_d, "d")


size = scaling[scaling.panel == "size"].sort_values("mb")
dens = scaling[scaling.panel == "density"].sort_values("genes_per_mb")


def max_ticks(ax, xs, ys, x_mid, y_mid):
    xmax, ymax = float(np.max(xs)), float(np.max(ys))
    ax.set_xticks([0, x_mid, xmax])
    ax.set_xticklabels(["0", f"{x_mid:g}", f"{xmax:.0f}"])
    ax.set_yticks([0, y_mid, ymax])
    ax.set_yticklabels(["0", f"{y_mid:g}", f"{ymax:.1f}"])


for ax, sub, xcol in ((ax_e, size, "mb"), (ax_f, dens, "genes_per_mb")):
    ax.plot(sub[xcol], sub.minutes, "-", color=CURVE, lw=0.9, alpha=0.35)
    curve(ax, sub[xcol], moving_average(sub.minutes), color=CURVE, lw=1.6)
    ax.set_xlim(0, sub[xcol].max() * 1.06)
    ax.set_ylim(0, sub.minutes.max() * 1.15)

max_ticks(ax_e, size.mb, size.minutes, 150, 2)
ax_e.set_xlabel("Genome size (Mb)"); ax_e.set_ylabel("Wall clock (min)")
key(ax_e, "e")
max_ticks(ax_f, dens.genes_per_mb, dens.minutes, 70, 0.5)
ax_f.set_xlabel("Genes per Mb")
key(ax_f, "f")


fig.canvas.draw()
pa, pb = ax_a.get_position(), ax_b.get_position()
gap = min(pa.y0, pb.y0) - max(ax_c.get_position().y1, ax_d.get_position().y1)
fig.text((pa.x0 + pb.x1) / 2, min(pa.y0, pb.y0) - 0.30 * gap,
         "Injected $K_s$", ha="center", va="top", fontsize=11)

handles, labels = ax_a.get_legend_handles_labels()
leg = ax_key.legend(handles, labels, loc="center left", ncol=1,
                    title="Heterozygosity\nof the assembly",
                    handletextpad=0.5, borderpad=0.0, handlelength=1.3,
                    borderaxespad=0.0, labelspacing=0.55,
                    bbox_to_anchor=(0.0, 0.5))
leg.get_title().set_fontsize(8.5)
leg.get_title().set_color(MUTED)

fig.savefig(snakemake.output.png)
fig.savefig(snakemake.output.pdf)
plt.close(fig)
