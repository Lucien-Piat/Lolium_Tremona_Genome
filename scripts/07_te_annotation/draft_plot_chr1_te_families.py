#!/usr/bin/env python3
# scripts/07_te_annotation/draft_plot_chr1_te_families.py  (plotting sif)
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BED="results/te_hite/gene_te/te.bed"; FAI="reference_data/lmultiflorum.tremona.fa.fai"
CHR="chr1"; WIN=1_000_000
OUT="results/te_hite/figures/chr1_te_families.pdf"
CENTRO_FAMILY=None          # set to None to disable the centromere band
os.makedirs(os.path.dirname(OUT), exist_ok=True)
None
SELECT=["chr5__LTR_612","chr6__LTR_648","chr2__LTR_970","chr3__LTR_774", "chr4__LTR_719",
        "chr4__TIR_410","chr4__TIR_212","chr6__TIR_168", "chr1__Helitron_1",
        "chr5__Denovo_Non_LTR_16","chr1__LTR_5"]
NAMES={}   # optional: {"chr6__LTR_648":"Centaurea", ...} to override auto labels

def category(c):
    if c=="LTR/Gypsy": return "Gypsy"
    if c=="LTR/Copia": return "Copia"
    if c.startswith("LTR/"): return "other LTR"
    if c.startswith("LINE"): return "LINE"
    if c.startswith("SINE"): return "SINE"
    if c.startswith("DNA/") or c.startswith("RC/"): return "DNA transposon"
    if c=="Unknown": return "Unknown"
    return "Other"
def code(c):
    m={"LTR/Gypsy":"RLG","LTR/Copia":"RLC","LTR/ERV":"RLE","LTR/Pao":"RLB",
       "LINE/L1":"RIL","LINE/RTE-RTE":"RIT","SINE/tRNA":"RST","DNA/PIF-Harbinger":"DTH",
       "DNA/CMC-EnSpm":"DTC","DNA/MULE":"DTM","DNA/hAT":"DTA","DNA/TcMar":"DTT",
       "RC/Helitron":"DHH","Unknown":"RLX"}
    return m.get(c,"XXX")
CCOL={"Gypsy":"#c0392b","Copia":"#27ae60","other LTR":"#16a085","LINE":"#8e44ad",
      "SINE":"#d2b4de","DNA transposon":"#2980b9","Unknown":"#95a5a6","Other":"#bdc3c7"}

df=pd.read_csv(BED, sep="\t", names=["chrom","start","end","family","klass","div","strand"])
df=df[(df.chrom==CHR)&(df.family.isin(SELECT))].copy()
missing=[f for f in SELECT if f not in set(df.family)]
if missing: print("WARNING: not on", CHR, ":", missing)
sizes={l.split('\t')[0]:int(l.split('\t')[1]) for l in open(FAI)}
clen=sizes[CHR]; nwin=clen//WIN+1
df["win"]=(df.start//WIN).clip(0,nwin-1).astype(int); df["len"]=df.end-df.start

mat,famcat,famcode,tot={},{},{},{}
for fam,g in df.groupby("family"):
    v=np.zeros(nwin); np.add.at(v, g.win.values, g.len.values)
    k=g.klass.mode().iat[0] if not g.klass.mode().empty else g.klass.iat[0]
    mat[fam]=v; tot[fam]=v.sum(); famcat[fam]=category(k); famcode[fam]=code(k)

# auto labels: code + rank within code
seen={}; auto={}
for fam in SELECT:
    if fam not in mat: continue
    c=famcode[fam]; seen[c]=seen.get(c,0)+1
    auto[fam]=NAMES.get(fam, f"{c}-{seen[c]}")

cen=None
if CENTRO_FAMILY in mat:
    v=mat[CENTRO_FAMILY]; w=np.where(v>=0.5*v.max())[0]
    if len(w): cen=(w.min()*WIN/1e6, (w.max()+1)*WIN/1e6)

order=[f for f in SELECT if f in mat]
x=np.arange(nwin)*WIN/1e6
fig,axes=plt.subplots(len(order),1,figsize=(10,0.66*len(order)),sharex=True)
axes=np.atleast_1d(axes)
for ax,fam in zip(axes,order):
    if cen: ax.axvspan(cen[0],cen[1],color="0.88",zorder=0)
    ax.fill_between(x, 0, mat[fam]/1e3, color=CCOL[famcat[fam]], linewidth=0)
    ax.set_yticks([]); ax.set_ylabel(f"{auto[fam]}\n{tot[fam]/1e6:.1f} Mb",
                                     rotation=0, ha="right", va="center", fontsize=7)
    for s in ("top","right","left"): ax.spines[s].set_visible(False)
axes[-1].set_xlabel(f"{CHR} position (Mb)")
fig.suptitle(f"TE family distribution along {CHR}", y=1.0)
fig.tight_layout(); fig.savefig(OUT, bbox_inches="tight"); print("Figure ->", OUT)