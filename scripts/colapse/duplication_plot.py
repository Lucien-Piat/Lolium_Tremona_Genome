#!/usr/bin/env python3
"""
Génération de la figure principale (fig_main.pdf) format "Poster" (haute, texte large)
à partir du fichier block_master.tsv pré-calculé.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8-white")
matplotlib.rcParams.update({
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "font.size": 25,
    "axes.titlesize": 24,
    "axes.labelsize": 25,
    "xtick.labelsize": 25,
    "ytick.labelsize": 25,
    "legend.fontsize": 25,
    "figure.titlesize": 26
})

MASTER_TSV = "results/dupclass/block_master.tsv"
OUTDIR = "results/dupclass"

KS_SPLIT = 0.2
ARTI_DEPTH = 0.8
ARTI_SAMEPOS = 0.10
SHARED_LABEL = "Shared polymorphism fraction"

COL_COMPLETE, COL_DUP = "#2e7d32", "#ef6c00"
COL_HIGH = "#c62828"
CALL_PAL = {"artifact": COL_DUP, "paralog": COL_COMPLETE}


def scatter_by_call(ax, master, x, y, vline, hline):
    for call, col in CALL_PAL.items():
        s = master[master["final_call"] == call].dropna(subset=[x, y])
        ax.scatter(s[x], s[y], c=col, s=s["n_pairs"] * 6, alpha=0.85,
                   edgecolor="white", lw=0.5, label=call)
    ax.axvline(vline, color="gray", ls="--", lw=2)
    ax.axhline(hline, color=COL_HIGH, ls=":", lw=2)
    ax.text(vline - 0.006, 0.97, "recent", transform=ax.get_xaxis_transform(),
            ha="right", va="top", fontsize=16, color="gray", style="italic")
    ax.text(vline + 0.006, 0.97, "old", transform=ax.get_xaxis_transform(),
            ha="left", va="top", fontsize=16, color="gray", style="italic")

def leaf_counts(master):
    old_par = int((master["median_ks"] >= KS_SPLIT).sum())
    recent = master["median_ks"] < KS_SPLIT
    arti = int((master["final_call"] == "artifact").sum())
    recent_par = int((recent & (master["final_call"] == "paralog")).sum())
    return old_par, arti, recent_par

def draw_decision_tree(ax, counts):
    old_par, arti, recent_par = counts
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("C.", loc="left", fontweight="bold")
    def box(x, y, txt, color, leaf=False):
        ax.text(x, y, txt, ha="center", va="center", fontsize=18,
                color="white" if leaf else "black", fontweight="medium",
                bbox=dict(boxstyle="round,pad=0.8", fc=color if leaf else "white",
                          ec="white" if leaf else color, lw=2))
    box(0.5, 0.9, f"Median Ks \u2265 {KS_SPLIT}", "black")
    box(0.82, 0.55, "paralog", CALL_PAL["paralog"], True)
    box(0.30, 0.55, f"Depth < {ARTI_DEPTH}\nor SharedPoly > {ARTI_SAMEPOS}", "black")
    box(0.13, 0.15, "artifact", CALL_PAL["artifact"], True)
    box(0.50, 0.15, "paralog", CALL_PAL["paralog"], True)
    
    ax.text(0.82, 0.43, f"n={old_par}", ha="center", fontsize=22, color="black")
    ax.text(0.13, 0.03, f"n={arti}", ha="center", fontsize=22, color="black")
    ax.text(0.50, 0.03, f"n={recent_par}", ha="center", fontsize=22, color="black")
    
    arr = dict(arrowstyle="->", color="gray", lw=2)
    ax.annotate("", xy=(0.78, 0.6), xytext=(0.56, 0.86), arrowprops=arr); ax.text(0.71, 0.76, "yes", fontsize=22, color="black")
    ax.annotate("", xy=(0.33, 0.62), xytext=(0.45, 0.86), arrowprops=arr); ax.text(0.34, 0.76, "no", fontsize=22, color="black")
    ax.annotate("", xy=(0.15, 0.21), xytext=(0.26, 0.49), arrowprops=arr); ax.text(0.14, 0.37, "yes", fontsize=22, color="black")
    ax.annotate("", xy=(0.48, 0.21), xytext=(0.34, 0.49), arrowprops=arr); ax.text(0.45, 0.37, "no", fontsize=22, color="black")

def draw_upset_split(fig, gs_slot, master):
    sig_names = [f"Ks < {KS_SPLIT}", f"Depth < {ARTI_DEPTH}", f"Shared > {ARTI_SAMEPOS}"]
    
    s_rec = master["median_ks"] < KS_SPLIT
    s_dep = master.get("depth_ratio_min", pd.Series(np.nan, index=master.index)) < ARTI_DEPTH
    s_pol = master.get("samepos_frac", pd.Series(np.nan, index=master.index)) > ARTI_SAMEPOS

    df = pd.DataFrame({
        sig_names[0]: s_rec.fillna(False),
        sig_names[1]: s_dep.fillna(False),
        sig_names[2]: s_pol.fillna(False),
        "Call": master["final_call"]
    })

    inner_main = gs_slot.subgridspec(1, 2, wspace=0.4)
    
    grps, global_max_y, max_cols = {}, 0, 0
    
    for call in ["artifact", "paralog"]:
        sub_df = df[df["Call"] == call]
        grp = sub_df.groupby(sig_names).size().reset_index(name="count")
        grp = grp[grp["count"] > 0].sort_values("count", ascending=False)
        grps[call] = grp
        if len(grp) > 0:
            global_max_y = max(global_max_y, grp["count"].max())
            max_cols = max(max_cols, len(grp))

    for idx, call in enumerate(["artifact", "paralog"]):
        grp = grps[call]
        sub_gs = inner_main[idx].subgridspec(2, 1, height_ratios=[2.5, 1], hspace=0.1)
        ax_bar = fig.add_subplot(sub_gs[0])
        ax_mat = fig.add_subplot(sub_gs[1], sharex=ax_bar)

        if len(grp) == 0:
            ax_bar.axis("off"); ax_mat.axis("off")
            continue

        x = np.arange(len(grp))
        y = grp["count"].values
        color = CALL_PAL.get(call, "gray")
        ax_bar.bar(x, y, color=color, edgecolor="white", lw=1, width=0.6)

        for i, val in enumerate(y):
            ax_bar.text(i, val + global_max_y * 0.02, str(int(val)), 
                        ha='center', va='bottom', fontsize=16)

        ax_bar.set_title(f"D1." if idx == 0 else f"D2.", loc="left", fontweight="bold", fontsize=24)
        
        for sp in ["top", "right", "bottom", "left"]: 
            ax_bar.spines[sp].set_visible(False)
        
        ax_bar.set_xticks([])
        ax_bar.set_ylim(0, global_max_y * 1.15)
        ax_bar.set_xlim(-0.5, max_cols - 0.5)
        
        if idx == 0: ax_bar.set_ylabel("Blocks")
        
        ax_mat.set_ylim(-0.5, 2.5)
        ax_mat.set_yticks([0, 1, 2])
        ax_mat.set_yticklabels(sig_names[::-1], fontsize=25)
        
        for spine in ax_mat.spines.values(): spine.set_visible(False)
        ax_mat.set_xticks([])
        ax_mat.tick_params(axis='y', length=0)

        for i, (_, row) in enumerate(grp.iterrows()):
            active_y = []
            for j, sig in enumerate(sig_names):
                y_pos = 2 - j 
                if row[sig]:
                    active_y.append(y_pos)
                    ax_mat.plot(i, y_pos, marker='o', color='#333333', markersize=14)
                else:
                    ax_mat.plot(i, y_pos, marker='o', color='#e0e0e0', markersize=14)
            
            if len(active_y) > 1:
                ax_mat.plot([i, i], [min(active_y), max(active_y)], color='#333333', lw=4)

        for y_pos in [0, 1, 2]:
            ax_mat.axhline(y_pos, color='#e0e0e0', lw=2, zorder=-1)

def fig_main(master, outdir):
    fig = plt.figure(figsize=(12, 22), constrained_layout=True)
    
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.8])

    ax = fig.add_subplot(gs[0, 0])
    scatter_by_call(ax, master, "median_ks", "depth_ratio_min", KS_SPLIT, ARTI_DEPTH)
    ax.set_xlabel("Median Ks"); ax.set_ylabel("Min depth / median")
    ax.set_title("A.", loc="left", fontweight="bold"); ax.legend(frameon=False)

    ax = fig.add_subplot(gs[0, 1])
    scatter_by_call(ax, master, "median_ks", "samepos_frac", KS_SPLIT, ARTI_SAMEPOS)
    ax.set_xlabel("Median Ks"); ax.set_ylabel(SHARED_LABEL)
    ax.set_title("B.", loc="left", fontweight="bold")

    draw_decision_tree(fig.add_subplot(gs[1, :]), leaf_counts(master))
    
    draw_upset_split(fig, gs[2, :], master)

    outpath = f"{outdir}/fig_main_poster.pdf"
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"Graphique généré avec succès : {outpath}")

if __name__ == "__main__":
    import os
    if not os.path.exists(MASTER_TSV):
        print(f"Erreur : Le fichier {MASTER_TSV} est introuvable. Exécutez le script depuis le bon dossier parent.")
    else:
        print(f"Chargement des données depuis {MASTER_TSV}...")
        df = pd.read_csv(MASTER_TSV, sep="\t")
        fig_main(df, OUTDIR)