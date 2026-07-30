"""
Synteny block Ks distributions converted to divergence time (Mya),
one stacked subplot per species.
T = Ks / (2 * mu)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# =====================================================================
# CONTROLS
# =====================================================================
X_MIN          = 0
X_MAX          = 200        # Mya
BIN_WIDTH      = 2.5        # Mya per bin

FIG_WIDTH      = 13
FIG_HEIGHT     = 14

COLOR_INTER    = "#4FA3DC"
COLOR_INTRA    = "#C846A8"
EDGE_COLOR     = "white"
EDGE_WIDTH     = 0.6

BASE_DIR       = "results/02_synteny"
FILE_NAME      = "blocks_ks.tsv"
OUTPUT         = "synteny_time_histograms.pdf"

# (mu in sub/site/year, citation). All from De La Torre et al. 2017, Table 1.
MUTATION_RATES = {
    "tremona":      (5.76e-9, "De La Torre et al. 2017"),
    "rabiosa":      (5.76e-9, "De La Torre et al. 2017"),
    "sikem":        (5.76e-9, "De La Torre et al. 2017"),
    "paraquat":     (5.76e-9, "De La Torre et al. 2017"),
    "perenne":      (5.76e-9, "De La Torre et al. 2017"),
    "brachypodium": (5.76e-9, "De La Torre et al. 2017"),
    "oryza":        (5.76e-9, "De La Torre et al. 2017"),
    "ananas":       (5.76e-9, "De La Torre et al. 2017"),
    "arabido":      (1.50e-8, "Koch et al. 2000"),
    "amborella":    (2.00e-9, "approx., slow basal angiosperm"),
}

# top to bottom order
SPECIES = [
    ("tremona",      "Lolium multiflorum Tremona"),
    ("paraquat",     "Lolium multiflorum Brunharo"),
    ("rabiosa",      "Lolium multiflorum Rabiosa"),
    ("sikem",        "Lolium multiflorum Sikem"),
    ("perenne",      "Lolium perenne Kyuss"),
    ("brachypodium", "Brachipodium dystachion"),
    ("oryza",        "Oryza sativa"),
    ("ananas",       "Ananas comosus"),
    ("arabido",      "Arabidopsis thaliana"),
    ("amborella",    "Amborella trichopoda"),
]

# =====================================================================
# LOAD + SHARED Y MAX
# =====================================================================
bins       = np.arange(X_MIN, X_MAX + BIN_WIDTH, BIN_WIDTH)
data_cache = {}
y_max      = 0

for folder, _ in SPECIES:
    path = os.path.join(BASE_DIR, folder, FILE_NAME)
    if not os.path.exists(path) or folder not in MUTATION_RATES:
        continue
    mu = MUTATION_RATES[folder][0]
    df = pd.read_csv(path, sep="\t")
    intra_mask = df["q_chr"] == df["t_chr"]
    inter_t = df.loc[~intra_mask, "median_ks"].to_numpy() / (2.0 * mu * 1e6)
    intra_t = df.loc[ intra_mask, "median_ks"].to_numpy() / (2.0 * mu * 1e6)
    h_inter, _ = np.histogram(inter_t, bins=bins)
    h_intra, _ = np.histogram(intra_t, bins=bins)
    stacked    = h_inter + h_intra
    if len(stacked):
        y_max = max(y_max, int(stacked.max()))
    data_cache[folder] = (inter_t, intra_t)

Y_MAX = int(np.ceil(y_max * 1.15)) if y_max else 10

# =====================================================================
# FIGURE
# =====================================================================
fig, axes = plt.subplots(
    nrows=len(SPECIES), ncols=1,
    figsize=(FIG_WIDTH, FIG_HEIGHT),
    sharex=True, sharey=True,
)

for ax, (folder, label) in zip(axes, SPECIES):
    if folder not in data_cache:
        ax.text(0.5, 0.5, f"missing: results/02_synteny/{folder}/{FILE_NAME}",
                transform=ax.transAxes, ha="center", va="center",
                color="grey", fontsize=9)
    else:
        inter_t, intra_t = data_cache[folder]
        mu, citation     = MUTATION_RATES[folder]

        ax.hist([inter_t, intra_t], bins=bins, stacked=True,
                color=[COLOR_INTER, COLOR_INTRA],
                edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH)

        ax.text(0.012, 0.92, label, transform=ax.transAxes,
                fontsize=20, fontweight="bold", va="top")

        mu_str = f"$\\mu$ = {mu:.2e}"
        handles = [
            Patch(facecolor=COLOR_INTER, edgecolor=EDGE_COLOR,
                  linewidth=EDGE_WIDTH,
                  label=f"Inter (n = {len(inter_t)})"),
            Patch(facecolor=COLOR_INTRA, edgecolor=EDGE_COLOR,
                  linewidth=EDGE_WIDTH,
                  label=f"Intra (n = {len(intra_t)})"),
            Patch(facecolor="none", edgecolor="none", label=mu_str),
        ]
        ax.legend(handles=handles, loc="upper right",
                  frameon=False, fontsize=13,
                  handlelength=1.2, labelspacing=0.5)

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(0, 20)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)

axes[-1].set_xlabel("Divergence time (Mya),  T = Ks / (2 $\\mu$)", fontsize=20)
fig.text(0.015, 0.5, "Number of synteny blocks",
         va="center", rotation="vertical", fontsize=20)

plt.tight_layout(rect=[0.04, 0, 1, 1])
plt.savefig(OUTPUT, format="pdf", bbox_inches="tight")
plt.show()