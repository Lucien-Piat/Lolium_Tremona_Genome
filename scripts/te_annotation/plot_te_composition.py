#!/usr/bin/env python3
# scripts/te_annotation/plot_te_composition.py  (plotting sif)
import os
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLASSTAB="results/te_hite/tremona_TE.class_table.tsv"
PART="results/te_hite/genome_partition.tsv"
OUT="results/te_hite/figures/te_composition_bar.pdf"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

PAL={"LTR/Gypsy":"#c0392b","LTR/Copia":"#27ae60","LTR/ERV":"#16a085","LTR/Pao":"#1abc9c",
     "LINE/L1":"#8e44ad","LINE/RTE-RTE":"#c39bd3","SINE/tRNA":"#d2b4de",
     "DNA/PIF-Harbinger":"#2980b9","DNA/CMC-EnSpm":"#3498db","DNA/MULE":"#5dade2",
     "DNA/hAT":"#85c1e9","DNA/TcMar":"#2e4053","RC/Helitron":"#aed6f1"}

d=pd.read_csv(CLASSTAB, sep="\t"); d.columns=["klass","frag","bp","pct"]
d=d.sort_values("bp", ascending=False)
p=pd.read_csv(PART, sep="\t"); part=dict(zip(p["component"], p["bp"]))
genome=float(part["genome_total"]); coding=float(part["coding_nonTE"])

segs=[]; minor=0.0
for _,r in d.iterrows():
    if r.pct>=0.2: segs.append((r.klass, r.bp, PAL.get(r.klass,"#999999")))
    else: minor+=r.bp
if minor>0: segs.append(("Other TE classes", minor, "#cfd2d4"))
te_total=d.bp.sum()
segs.append(("Unannotated", max(genome-te_total-coding,0), "#ffffff"))   # before coding
segs.append(("Coding DNA", coding, "#f1c40f"))                            # last, separated

fig, ax=plt.subplots(figsize=(11,3.2))
left=0.0
for lab,bp,col in segs:
    g=bp/1e9
    ax.barh(0, g, left=left, color=col, edgecolor="#bbbbbb", linewidth=0.4)
    if g>=0.06: ax.text(left+g/2, 0, f"{g:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if col not in ("#ffffff","#f1c40f","#cfd2d4") else "#333")
    left+=g
ax.set_yticks([0]); ax.set_yticklabels(["Tremona"]); ax.set_ylim(-0.6,0.6)
ax.set_xlabel("Genome (Gb)"); ax.set_title("Genome composition")
handles=[plt.Rectangle((0,0),1,1,color=c,ec="#bbbbbb") for _,_,c in segs]
labels=[f"{l} ({bp/genome*100:.1f}%)" for l,bp,_ in segs]
ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5,-0.85),
          ncol=5, fontsize=8, frameon=False, handlelength=1, columnspacing=1)
fig.subplots_adjust(bottom=0.45)
fig.savefig(OUT, bbox_inches="tight"); print("Figure ->", OUT)