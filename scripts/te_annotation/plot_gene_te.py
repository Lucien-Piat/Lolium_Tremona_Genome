#!/usr/bin/env python3
# scripts/te_annotation/plot_gene_te.py  (plotting sif, after 09b)
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

DIST="results/te_hite/gene_te/te_gene_distance.tsv"
OUT="results/te_hite/figures/gene_te.pdf"; NEAR=2000; ALPHA=0.85
os.makedirs(os.path.dirname(OUT), exist_ok=True)

def category(c):
    if c=="LTR/Gypsy": return "Gypsy"
    if c=="LTR/Copia": return "Copia"
    if c.startswith("LTR/"): return "other LTR"
    if c.startswith("LINE"): return "LINE"
    if c.startswith("SINE"): return "SINE"
    if c.startswith("DNA/") or c.startswith("RC/"): return "DNA transposon"
    if c=="Unknown": return "Unknown"
    return "Other"
CATS=["Gypsy","Copia","other LTR","LINE","DNA transposon","SINE","Unknown","Other"]
CCOL={"Gypsy":"#c0392b","Copia":"#27ae60","other LTR":"#16a085","LINE":"#8e44ad",
      "DNA transposon":"#2980b9","SINE":"#d2b4de","Unknown":"#95a5a6","Other":"#bdc3c7"}
FACE={c: to_rgba(CCOL[c], ALPHA) for c in CATS}

df=pd.read_csv(DIST, sep="\t", names=["distance","klass","length"])
df["cat"]=df.klass.map(category)

# left: stacked density vs distance, log-spaced bins (in gene -> first bin)
edges=[0,100,200,500,1000,2000,5000,10000,20000,30000, np.inf]
cen=[]
for i in range(len(edges)-1):
    lo,hi=edges[i],edges[i+1]
    cen.append((50 if lo==0 else (40000 if hi==np.inf else (lo*hi)**0.5))/1000.0)
df["bin"]=pd.cut(df.distance, bins=edges, right=False, labels=False)
g=df.groupby(["bin","cat"])["length"].sum().unstack(fill_value=0).reindex(range(len(cen)), fill_value=0)
frac=g.div(g.sum(axis=1), axis=0).fillna(0)
ys=[frac[c].values if c in frac.columns else np.zeros(len(cen)) for c in CATS]

# right: enrichment within NEAR vs genome-wide share
allbp=df.groupby("cat")["length"].sum()
nearbp=df[df.distance<=NEAR].groupby("cat")["length"].sum().reindex(allbp.index).fillna(0)
enr=np.log2(((nearbp/nearbp.sum())/(allbp/allbp.sum())).replace(0,np.nan)).dropna().sort_values()

fig,axes=plt.subplots(1,2,figsize=(12,4.6))
axes[0].stackplot(cen, ys, colors=[FACE[c] for c in CATS], edgecolor="black", linewidth=0.2)
axes[0].set_xscale("log"); axes[0].set_xlim(0.04,45); axes[0].set_ylim(0,1)
axes[0].set_xticks([0.1,0.5,1,5,10,30]); axes[0].set_xticklabels(["0.1","0.5","1","5","10","30"])
axes[0].set_xlabel("distance to nearest gene (kb, log scale)")
axes[0].set_ylabel("fraction of TE bp"); axes[0].set_title("TE composition vs distance to genes")

cols=[to_rgba(CCOL.get(c,"#999999"), ALPHA) for c in enr.index]
axes[1].barh(range(len(enr)), enr.values, color=cols, edgecolor="black", linewidth=0.2)
axes[1].axvline(0, color="k", lw=0.8)
axes[1].set_yticks(range(len(enr))); axes[1].set_yticklabels(enr.index, fontsize=8)
axes[1].set_xlabel("log2 enrichment within 2 kb of a gene")
axes[1].set_title("Gene-proximal enrichment by class")
for s in ("top","right"): axes[1].spines[s].set_visible(False)

handles=[plt.Rectangle((0,0),1,1, facecolor=FACE[c], edgecolor="black", linewidth=0.2) for c in CATS]
axes[1].legend(handles, CATS, loc="lower right", fontsize=7, frameon=False, ncol=2)
fig.tight_layout(); fig.savefig(OUT, bbox_inches="tight"); print("Figure ->", OUT)