#!/usr/bin/env python3
"""
Classification of Tremona self-synteny blocks (collapse artifact vs genuine
paralog) and masking of one copy per artifact block.
"""

import argparse
import re
import subprocess
from pathlib import Path

import numpy as np # type: ignore
import pandas as pd # type: ignore
import pysam # type: ignore
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore
from matplotlib.patches import Patch, Circle, Rectangle # type: ignore
import seaborn as sns # type: ignore

plt.style.use("seaborn-v0_8-white")
matplotlib.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": False})

CONFIG = {
    "vcf":   "reference_data/cohort.snps.vcf.gz",
    "hwe":   "collapse_diag/trem.hwe",
    "bam":   "mapping/TREM1.dedup.bam",
    "ks":    "results/02_synteny/tremona/blocks_ks.tsv",
    "collinearity": "results/02_synteny/tremona/tremona.collinearity",
    "gff":   "results/02_synteny/tremona/tremona.gff",
    "gff3":  "results/02_synteny/tremona/annotation.gff3",
    "axt":   "results/02_synteny/tremona/kaks_dapnv60w/pairs.axt",
    "fai":   "reference_data/lmultiflorum.tremona.fa.fai",
    "busco": "reference_data/lmultiflorum.tremona_full_table_busco_format.tsv",
    "presence_matrix": "results/02_duplication_sharing/presence_matrix.tsv",
    "outdir": "results/03_duplication_classification",
}
MOSDEPTH = "mosdepth"
MAPQ_HI = 20
WIN = 50000
KS_SPLIT = 0.2
ARTI_DEPTH = 0.8
ARTI_SAMEPOS = 0.10
MIN_EITHER = 15
SMOOTH = 100
COLLAPSE_DEPTH = 1.5
CHROMS = [f"chr{i}" for i in range(1, 8)]
SHARED_LABEL = "Shared polymorphism fraction"


COL_INTER, COL_INTRA = "#1976d2", "#9c27b0"
COL_COMPLETE, COL_DUP, COL_FRAG, COL_MISSING = "#2e7d32", "#ef6c00", "#fbc02d", "#9e9e9e"
COL_HIGH, COL_LOW = "#c62828", "#1565c0"   # above / below baseline
COL_MASK = "#2b2d42"
LOC_PAL = {"intra": COL_INTRA, "inter": COL_INTER}
BUSCO_PAL = {"S": COL_COMPLETE, "D": COL_DUP, "F": COL_FRAG, "M": COL_MISSING}
CALL_PAL = {"artifact": COL_DUP, "paralog": COL_COMPLETE}

def norm_chr(c):
    c = str(c)
    return "chr" + c[2:] if c.startswith("tr") else c

def log(msg):
    print(msg, flush=True)

def prog(label, i, n, every=20):
    if i % every == 0 or i == n:
        log(f"    {label} {i}/{n}")

def cached(path, fn, force=False):
    path = Path(path)
    if path.exists() and not force:
        log(f"  [cache] {path.name}")
        return pd.read_csv(path, sep="\t")
    df = fn()
    df.to_csv(path, sep="\t", index=False)
    return df

# Data loading

def parse_collinearity(path):
    blocks, cur = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("## Alignment"):
                cur = int(line.split("Alignment")[1].split(":")[0]); blocks[cur] = ([], [])
            elif line.startswith("#") or cur is None:
                continue
            else:
                p = line.strip().split("\t")
                if len(p) >= 3:
                    blocks[cur][0].append(p[1]); blocks[cur][1].append(p[2])
    return blocks

def load_genes():
    g = pd.read_csv(CONFIG["gff"], sep="\t", header=None, names=["chr", "gene", "start", "end"])
    g["chr"] = g["chr"].map(norm_chr)
    g["start"], g["end"] = g[["start", "end"]].min(axis=1), g[["start", "end"]].max(axis=1)
    return g

def genes_in(genes, chrom, start, end):
    return genes[(genes["chr"] == chrom) & (genes["end"] >= start) & (genes["start"] <= end)]

def load_blocks(genes):
    ks = pd.read_csv(CONFIG["ks"], sep="\t")
    for c in ("q_chr", "t_chr"):
        ks[c] = ks[c].map(norm_chr)
    coll = parse_collinearity(CONFIG["collinearity"])
    gpos = {r["gene"]: (r["chr"], int(r["start"]), int(r["end"])) for _, r in genes.iterrows()}
    rows, missing = [], 0
    for _, r in ks.iterrows():
        bid = int(r["block_id"])
        if bid not in coll:
            missing += 1; continue
        qc = [gpos[g] for g in coll[bid][0] if g in gpos]
        tc = [gpos[g] for g in coll[bid][1] if g in gpos]
        if not qc or not tc:
            missing += 1; continue
        q_chr = pd.Series([c for c, _, _ in qc]).mode().iloc[0]
        t_chr = pd.Series([c for c, _, _ in tc]).mode().iloc[0]
        qs = [s for c, s, _ in qc if c == q_chr]; qe = [e for c, _, e in qc if c == q_chr]
        ts = [s for c, s, _ in tc if c == t_chr]; te = [e for c, _, e in tc if c == t_chr]
        rows.append({"block_id": bid, "q_chr": q_chr, "q_start": min(qs), "q_end": max(qe),
                     "t_chr": t_chr, "t_start": min(ts), "t_end": max(te),
                     "n_genes": len(coll[bid][0]), "n_pairs": r["n_pairs"], "median_ks": r["median_ks"]})
    df = pd.DataFrame(rows)
    if missing:
        log(f"  load_blocks : {missing} blocks without coordinates")
    df["is_intra"] = df["q_chr"] == df["t_chr"]
    df["loc"] = np.where(df["is_intra"], "intra", "inter")
    df["label_ks"] = np.where(df["median_ks"] < KS_SPLIT, "recent", "old")
    return df

def run_mosdepth(prefix, mapq, by, force=False):
    out = f"{prefix}.regions.bed.gz"
    if Path(out).exists() and not force:
        return out
    subprocess.run([MOSDEPTH, "-x", "-n", "-Q", str(mapq), "-t", "4", "--by", str(by),
                    prefix, CONFIG["bam"]], check=True)
    return out

def read_regions(path, named):
    cols = ["chr", "start", "end", "name", "depth"] if named else ["chr", "start", "end", "depth"]
    return pd.read_csv(path, sep="\t", header=None, names=cols)

def write_block_bed(df, path):
    rows = []
    for _, r in df.iterrows():
        rows.append((r["q_chr"], int(r["q_start"]), int(r["q_end"]), f"b{r['block_id']}:q"))
        rows.append((r["t_chr"], int(r["t_start"]), int(r["t_end"]), f"b{r['block_id']}:t"))
    pd.DataFrame(rows, columns=["chr", "start", "end", "name"]).to_csv(path, sep="\t", header=False, index=False)

def coverage_metrics(df, outdir, force=False):
    cov_path, win_path = outdir / "blocks_coverage.tsv", outdir / "genome_depth_win.tsv"
    if cov_path.exists() and win_path.exists() and not force:
        log("  [cache] blocks_coverage.tsv, genome_depth_win.tsv")
        return pd.read_csv(cov_path, sep="\t"), pd.read_csv(win_path, sep="\t")
    bed = outdir / "blocks.bed"; write_block_bed(df, bed)
    win = read_regions(run_mosdepth(str(outdir / "genome_q20"), MAPQ_HI, WIN, force), named=False)
    med = win["depth"].median()
    win = win[win["chr"].isin(CHROMS)].copy()
    win["win"], win["depth_ratio"] = win["start"], win["depth"] / med
    win = win[["chr", "win", "depth_ratio"]]
    q0 = read_regions(run_mosdepth(str(outdir / "blocks_q0"), 0, str(bed), force), named=True)
    q20 = read_regions(run_mosdepth(str(outdir / "blocks_q20"), MAPQ_HI, str(bed), force), named=True)
    m = q0.merge(q20, on=["chr", "start", "end", "name"], suffixes=("_q0", "_q20"))
    m["bid"] = m["name"].str.split(":").str[0].str[1:].astype(int)
    m["side"] = m["name"].str.split(":").str[1]
    m["depth_ratio"] = m["depth_q20"] / med
    m["multimap"] = 1 - (m["depth_q20"] / m["depth_q0"].replace(0, np.nan))
    qs = m[m["side"] == "q"].set_index("bid"); ts = m[m["side"] == "t"].set_index("bid")
    cov = pd.DataFrame(index=sorted(m["bid"].unique()))
    cov["depth_ratio_q"], cov["depth_ratio_t"] = qs["depth_ratio"], ts["depth_ratio"]
    cov["multimap_q"], cov["multimap_t"] = qs["multimap"], ts["multimap"]
    cov["depth_ratio_min"] = cov[["depth_ratio_q", "depth_ratio_t"]].min(axis=1)
    cov["depth_ratio_mean"] = cov[["depth_ratio_q", "depth_ratio_t"]].mean(axis=1)
    cov["multimap_mean"] = cov[["multimap_q", "multimap_t"]].mean(axis=1)
    cov = cov.reset_index().rename(columns={"index": "block_id"})
    cov.to_csv(cov_path, sep="\t", index=False); win.to_csv(win_path, sep="\t", index=False)
    return cov, win

def snp_in_genes(df, genes):
    vcf = pysam.VariantFile(CONFIG["vcf"]); out = []; n = len(df)
    for i, (_, r) in enumerate(df.iterrows(), 1):
        n_snp = n_het = glen = 0
        for chrom, s, e in [(r["q_chr"], r["q_start"], r["q_end"]), (r["t_chr"], r["t_start"], r["t_end"])]:
            for _, gene in genes_in(genes, chrom, int(s), int(e)).iterrows():
                glen += gene["end"] - gene["start"] + 1
                try:
                    it = vcf.fetch(chrom, int(gene["start"]), int(gene["end"]))
                except ValueError:
                    continue
                for rec in it:
                    if len(rec.ref) != 1 or any(len(a) != 1 for a in rec.alts or []):
                        continue
                    n_snp += 1
                    for smp in rec.samples.values():
                        gt = smp.get("GT")
                        if gt and None not in gt and gt[0] != gt[1]:
                            n_het += 1
        out.append((r["block_id"], n_snp, n_het, glen, n_snp / glen * 1000 if glen else np.nan))
        prog("snp_in_genes", i, n)
    return pd.DataFrame(out, columns=["block_id", "n_snp_genic", "n_het_genic", "genic_bp", "snp_per_kb"])

def load_hwe():
    h = pd.read_csv(CONFIG["hwe"], sep="\t")
    h.columns = [c.strip() for c in h.columns]
    h["chr"] = h[h.columns[0]].map(norm_chr)
    h["pos"] = h[h.columns[1]].astype(int)
    h["ho"] = h[h.columns[2]].map(lambda s: float(str(s).split("/")[1]))
    h["he"] = h[h.columns[3]].map(lambda s: float(str(s).split("/")[1]))
    return h[["chr", "pos", "ho", "he"]]

def fis_per_block(df, hwe):
    out = []
    for _, r in df.iterrows():
        mask = (((hwe["chr"] == r["q_chr"]) & (hwe["pos"].between(r["q_start"], r["q_end"])))
                | ((hwe["chr"] == r["t_chr"]) & (hwe["pos"].between(r["t_start"], r["t_end"]))))
        sub = hwe[mask]; he = sub["he"].sum()
        out.append((r["block_id"], 1 - sub["ho"].sum() / he if he > 0 else np.nan, len(sub)))
    return pd.DataFrame(out, columns=["block_id", "fis_block", "n_sites_hwe"])

def fis_windows(hwe):
    h = hwe.copy(); h["win"] = (h["pos"] // WIN) * WIN
    g = h.groupby(["chr", "win"]).agg(ho=("ho", "sum"), he=("he", "sum")).reset_index()
    g = g[g["he"] > 0].copy(); g["fis"] = 1 - g["ho"] / g["he"]
    return g[["chr", "win", "fis"]]

AXT_HEADER = re.compile(r"^B(\d+)__(.+?)__(.+)$")

def parse_axt(path):
    lines = [l.strip() for l in open(path) if l.strip()]
    recs, i = 0, 0
    recs_list = []
    while i < len(lines):
        m = AXT_HEADER.match(lines[i])
        if m and i + 2 < len(lines):
            recs_list.append((int(m.group(1)), m.group(2), m.group(3), lines[i + 1], lines[i + 2])); i += 3
        else:
            i += 1
    return recs_list

def parse_cds(gff3, needed):
    cds, gene2mrna = {}, {}
    with open(gff3) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attr = f[:9]
            ad = dict(kv.split("=", 1) for kv in attr.split(";") if "=" in kv)
            if ftype == "mRNA" and ad.get("Parent") and ad.get("ID"):
                gene2mrna.setdefault(ad["Parent"], ad["ID"])
            elif ftype == "CDS" and ad.get("Parent"):
                cds.setdefault(ad["Parent"], []).append((int(start), int(end), strand, chrom))

    def positions(segs):
        strand, chrom = segs[0][2], segs[0][3]
        segs = sorted(segs, key=lambda x: x[0], reverse=(strand == "-"))
        pos = []
        for s, e, _, _ in segs:
            pos.extend(range(e, s - 1, -1) if strand == "-" else range(s, e + 1))
        return chrom, pos

    out = {}
    for nid in needed:
        key = nid
        if key not in cds:
            stripped = re.sub(r"_\d+$", "", nid)
            if stripped in cds:
                key = stripped
            elif nid in gene2mrna and gene2mrna[nid] in cds:
                key = gene2mrna[nid]
            else:
                continue
        out[nid] = positions(cds[key])
    return out

def snpset_for(vcf, chrom, start, end, cache):
    k = (chrom, start, end)
    if k in cache:
        return cache[k]
    s = set()
    try:
        for rec in vcf.fetch(chrom, max(0, start - 1), end):
            if len(rec.ref) == 1 and rec.alts and all(len(a) == 1 for a in rec.alts):
                s.add(rec.pos)
    except ValueError:
        pass
    cache[k] = s
    return s

def same_position_metric(df):
    recs = parse_axt(CONFIG["axt"])
    cdsmap = parse_cds(CONFIG["gff3"], {g for _, g1, g2, _, _ in recs for g in (g1, g2)})
    vcf = pysam.VariantFile(CONFIG["vcf"])
    cache, per_block, skipped, n = {}, {}, 0, len(recs)
    for idx, (bid, g1, g2, s1, s2) in enumerate(recs, 1):
        prog("shared_polymorphism", idx, n, every=500)
        if g1 not in cdsmap or g2 not in cdsmap:
            skipped += 1; continue
        (c1, p1), (c2, p2) = cdsmap[g1], cdsmap[g2]
        u1, u2 = s1.replace("-", ""), s2.replace("-", "")
        if len(u1) != len(p1):
            if len(u1) < len(p1): p1 = p1[:len(u1)]
            else: skipped += 1; continue
        if len(u2) != len(p2):
            if len(u2) < len(p2): p2 = p2[:len(u2)]
            else: skipped += 1; continue
        snp1 = snpset_for(vcf, c1, min(p1), max(p1), cache)
        snp2 = snpset_for(vcf, c2, min(p2), max(p2), cache)
        i1 = i2 = both = either = 0
        for a, b in zip(s1, s2):
            if a != "-" and b != "-":
                p, q = p1[i1] in snp1, p2[i2] in snp2
                if p or q:
                    either += 1
                    if p and q: both += 1
            if a != "-": i1 += 1
            if b != "-": i2 += 1
        d = per_block.setdefault(bid, [0, 0, 0]); d[0] += both; d[1] += either; d[2] += 1
    if skipped:
        log(f"    shared_polymorphism : {skipped} pairs skipped")
    return pd.DataFrame([(b, v[0], v[1], v[2], v[0] / v[1] if v[1] else np.nan) for b, v in per_block.items()],
                        columns=["block_id", "samepos_both", "samepos_either", "samepos_pairs_used", "samepos_frac"])

def load_sharing():
    p = Path(CONFIG["presence_matrix"])
    if not p.exists():
        return None
    pm = pd.read_csv(p, sep="\t")
    if "block" not in pm.columns or "category" not in pm.columns:
        return None
    out = pm[["block", "category"]].rename(columns={"block": "block_id"})
    out["block_id"] = out["block_id"].astype(int)
    return out

# Decision, masking, BUSCO

def _final_call(r):
    if pd.notna(r["median_ks"]) and r["median_ks"] >= KS_SPLIT:
        return "paralog"
    depth, samepos = r.get("depth_ratio_min", np.nan), r.get("samepos_frac", np.nan)
    confirmed = (pd.notna(depth) and depth < ARTI_DEPTH) or (pd.notna(samepos) and samepos > ARTI_SAMEPOS)
    return "artifact" if confirmed else "paralog"

def leaf_counts(master):
    old_par = int((master["median_ks"] >= KS_SPLIT).sum())
    recent = master["median_ks"] < KS_SPLIT
    arti = int((master["final_call"] == "artifact").sum())
    recent_par = int((recent & (master["final_call"] == "paralog")).sum())
    return old_par, arti, recent_par

def masked_intervals(master, busco_loci=None):
    out = []
    for _, b in master[master["final_call"] == "artifact"].iterrows():
        qc, qs, qe = b["q_chr"], b["q_start"], b["q_end"]
        tc, ts, te = b["t_chr"], b["t_start"], b["t_end"]
        
        q_busco = 0
        t_busco = 0
        
        if busco_loci is not None:
            q_busco = len(busco_loci[(busco_loci["kind"] == "single") & 
                                     (busco_loci["chr"] == qc) & 
                                     (busco_loci["start"] <= qe) & 
                                     (busco_loci["end"] >= qs)])
            t_busco = len(busco_loci[(busco_loci["kind"] == "single") & 
                                     (busco_loci["chr"] == tc) & 
                                     (busco_loci["start"] <= te) & 
                                     (busco_loci["end"] >= ts)])

        # Prioritize masking the region with the LEAST complete (single) BUSCOs
        if q_busco < t_busco:
            ch, s, e = qc, qs, qe
        elif t_busco < q_busco:
            ch, s, e = tc, ts, te
        else:
            # Fallback to depth logic if they tie
            dq, dt = b.get("depth_ratio_q", np.nan), b.get("depth_ratio_t", np.nan)
            if pd.notna(dq) and pd.notna(dt) and dq < dt:
                ch, s, e = qc, qs, qe
            else:
                ch, s, e = tc, ts, te
                
        out.append((b["block_id"], ch, int(s), int(e)))
    return pd.DataFrame(out, columns=["block_id", "chr", "start", "end"])

def _busco_cols(bt):
    scol = next(c for c in bt.columns if "Status" in c)
    seqc = next(c for c in bt.columns if "Sequence" in c)
    sc = next(c for c in bt.columns if "Start" in c)
    ec = next(c for c in bt.columns if "End" in c)
    return scol, seqc, sc, ec

def load_busco_loci():
    p = Path(CONFIG["busco"])
    if not p.exists():
        return None
    bt = pd.read_csv(p, sep="\t")
    try:
        scol, seqc, sc, ec = _busco_cols(bt)
    except StopIteration:
        return None
    bt = bt.copy()
    bt["chr"] = bt[seqc].map(lambda x: norm_chr(x) if pd.notna(x) else x)
    bt["start"] = pd.to_numeric(bt[sc], errors="coerce")
    bt["end"] = pd.to_numeric(bt[ec], errors="coerce")
    st = bt[scol].astype(str)
    bt["kind"] = np.where(st.str.startswith("Dup"), "dup",
                 np.where(st.str.startswith("Complete"), "single", "other"))
    bt = bt[(bt["kind"].isin(["single", "dup"]))].dropna(subset=["start", "end"])
    return bt[bt["chr"].isin(CHROMS)][["chr", "start", "end", "kind"]]

def busco_recompute(masked):
    p = Path(CONFIG["busco"])
    if not p.exists():
        return None
    bt = pd.read_csv(p, sep="\t")
    try:
        scol, seqc, sc, ec = _busco_cols(bt)
    except StopIteration:
        return None
    idc = bt.columns[0]
    bt["_chr"] = bt[seqc].map(lambda x: norm_chr(x) if pd.notna(x) else x)
    bt[sc] = pd.to_numeric(bt[sc], errors="coerce")
    bt[ec] = pd.to_numeric(bt[ec], errors="coerce")
    st = bt[scol].astype(str)
    cond = [st.str.startswith("Complete"), st.str.startswith("Dup"), st.str.startswith("Frag")]
    bt["_b"] = np.select(cond, ["S", "D", "F"], default="M")
    status = bt.groupby(idc)["_b"].agg(lambda x: "D" if (x == "D").any() else x.iloc[0])
    tot = len(status)
    before = status.value_counts().to_dict()

    def in_masked(r):
        m = masked[(masked["chr"] == r["_chr"]) & (masked["start"] <= r[ec]) & (masked["end"] >= r[sc])]
        return len(m) > 0

    d2s = d2m = s2m = 0
    for did, stt in status.items():
        rows = bt[bt[idc] == did].dropna(subset=[sc, ec])
        if len(rows) == 0:
            continue
        remaining = len(rows) - sum(in_masked(r) for _, r in rows.iterrows())
        if stt == "D":
            if remaining <= 0: d2m += 1
            elif remaining == 1: d2s += 1
        elif stt == "S" and remaining <= 0:
            s2m += 1
    S0, D0, F0, M0 = (before.get(k, 0) for k in ("S", "D", "F", "M"))
    after = {"S": S0 + d2s - s2m, "D": D0 - d2s - d2m, "F": F0, "M": M0 + d2m + s2m}
    pct = lambda v: 100 * v / tot
    return {"tot": tot,
            "S_before": pct(S0), "D_before": pct(D0), "F_before": pct(F0), "M_before": pct(M0),
            "S_after": pct(after["S"]), "D_after": pct(after["D"]),
            "F_after": pct(after["F"]), "M_after": pct(after["M"])}

# drawing helpers
def scatter_by_call(ax, master, x, y, vline, hline):
    for call, col in CALL_PAL.items():
        s = master[master["final_call"] == call].dropna(subset=[x, y])
        ax.scatter(s[x], s[y], c=col, s=s["n_pairs"] * 3, alpha=0.85,
                   edgecolor="white", lw=0.5, label=call)
    ax.axvline(vline, color="gray", ls="--", lw=1)
    ax.axhline(hline, color=COL_HIGH, ls=":", lw=1)
    ax.text(vline - 0.006, 0.97, "recent", transform=ax.get_xaxis_transform(),
            ha="right", va="top", fontsize=8, color="gray", style="italic")
    ax.text(vline + 0.006, 0.97, "old", transform=ax.get_xaxis_transform(),
            ha="left", va="top", fontsize=8, color="gray", style="italic")

def draw_decision_tree(ax, counts):
    old_par, arti, recent_par = counts
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("C.", loc="left", fontweight="bold")
    def box(x, y, txt, color, leaf=False):
        ax.text(x, y, txt, ha="center", va="center", fontsize=9,
                color="white" if leaf else "black", fontweight="medium",
                bbox=dict(boxstyle="round,pad=0.5", fc=color if leaf else "white",
                          ec="white" if leaf else color, lw=1.5))
    box(0.5, 0.9, f"Median Ks \u2265 {KS_SPLIT}", "gray")
    box(0.82, 0.55, "paralog", CALL_PAL["paralog"], True)
    box(0.30, 0.55, f"Depth < {ARTI_DEPTH}\nor SharedPoly > {ARTI_SAMEPOS}", "gray")
    box(0.13, 0.15, "artifact", CALL_PAL["artifact"], True)
    box(0.50, 0.15, "paralog", CALL_PAL["paralog"], True)
    ax.text(0.82, 0.43, f"n={old_par}", ha="center", fontsize=8, color="gray")
    ax.text(0.13, 0.03, f"n={arti}", ha="center", fontsize=8, color="gray")
    ax.text(0.50, 0.03, f"n={recent_par}", ha="center", fontsize=8, color="gray")
    arr = dict(arrowstyle="->", color="gray", lw=1.2)
    ax.annotate("", xy=(0.78, 0.6), xytext=(0.56, 0.86), arrowprops=arr); ax.text(0.71, 0.76, "yes", fontsize=8, color="gray")
    ax.annotate("", xy=(0.33, 0.62), xytext=(0.45, 0.86), arrowprops=arr); ax.text(0.34, 0.76, "no", fontsize=8, color="gray")
    ax.annotate("", xy=(0.15, 0.21), xytext=(0.26, 0.49), arrowprops=arr); ax.text(0.14, 0.37, "yes", fontsize=8, color="gray")
    ax.annotate("", xy=(0.48, 0.21), xytext=(0.34, 0.49), arrowprops=arr); ax.text(0.45, 0.37, "no", fontsize=8, color="gray")

def draw_upset_split(fig, gs_slot, master):
    """Draws two side-by-side UpSet plots for artifacts vs paralogs using shared axes and threshold labels"""
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
    
    # global maximums
    grps = {}
    global_max_y = 0
    max_cols = 0
    
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
            ax_bar.axis("off")
            ax_mat.axis("off")
            continue

        x = np.arange(len(grp))
        y = grp["count"].values

        color = CALL_PAL.get(call, "gray")
        ax_bar.bar(x, y, color=color, edgecolor="white", lw=0.5, width=0.6)

        for i, val in enumerate(y):
            ax_bar.text(i, val + global_max_y * 0.02, str(int(val)), 
                        ha='center', va='bottom', fontsize=8)

        title = f"D1." if idx == 0 else f"D2."
        ax_bar.set_title(title, loc="left", fontweight="bold", fontsize=10)
        
        for sp in ["top", "right", "bottom", "left"]: 
            ax_bar.spines[sp].set_visible(False)
        
        ax_bar.set_xticks([])
        
        # SHARED AXES SCALING
        ax_bar.set_ylim(0, global_max_y * 1.15)
        ax_bar.set_xlim(-0.5, max_cols - 0.5)
        
        if idx == 0:
            ax_bar.set_ylabel("Blocks")
        
        ax_mat.set_ylim(-0.5, 2.5)
        ax_mat.set_yticks([0, 1, 2])
        ax_mat.set_yticklabels(sig_names[::-1], fontsize=8)
        
        for spine in ax_mat.spines.values(): 
            spine.set_visible(False)
            
        ax_mat.set_xticks([])
        ax_mat.tick_params(axis='y', length=0)

        for i, (_, row) in enumerate(grp.iterrows()):
            active_y = []
            for j, sig in enumerate(sig_names):
                y_pos = 2 - j 
                if row[sig]:
                    active_y.append(y_pos)
                    ax_mat.plot(i, y_pos, marker='o', color='#333333', markersize=6)
                else:
                    ax_mat.plot(i, y_pos, marker='o', color='#e0e0e0', markersize=6)
            
            if len(active_y) > 1:
                ax_mat.plot([i, i], [min(active_y), max(active_y)], color='#333333', lw=2)

        for y_pos in [0, 1, 2]:
            ax_mat.axhline(y_pos, color='#e0e0e0', lw=1, zorder=-1)
            
def draw_busco_bar(ax, br):
    ax.set_title("E.", loc="left", fontweight="bold")
    if not br:
        ax.axis("off"); return
    
    seg = [("S", "Complete"), ("D", "Duplicated"), ("F", "Fragmented"), ("M", "Missing")]
    rows = [("Before", "_before", 1), ("After", "_after", 0)]
    h = 0.5
    
    coords = {"Before": {}, "After": {}}
    
    for name, suf, yi in rows:
        left = 0.0
        for key, _ in seg:
            w = br[key + suf]
            ax.barh(yi, w, h, left=left, color=BUSCO_PAL[key], zorder=2)
            coords[name][key] = (left, left + w)
            
            if w >= 2.0:
                ax.text(left + w / 2, yi, f"{w:.1f}", va="center", ha="center",
                        color="white", fontsize=9, fontweight="bold", zorder=3)
            left += w

    y_bottom_before = 1 - h/2
    y_top_after = 0 + h/2
    
    # Sankey style 
    for key, _ in seg:
        x1_b, x2_b = coords["Before"][key]
        x1_a, x2_a = coords["After"][key]
        
        poly = plt.Polygon(
            [[x1_b, y_bottom_before], [x2_b, y_bottom_before], [x2_a, y_top_after], [x1_a, y_top_after]],
            facecolor=BUSCO_PAL[key], alpha=0.15, zorder=1
        )
        ax.add_patch(poly)

    xb_mid = sum(coords["Before"]["S"]) / 2
    xa_mid = sum(coords["After"]["S"]) / 2
    
    diff = br['S_after'] - br['S_before']
    sign = "+" if diff >= 0 else ""
    
    ax.text((xb_mid + xa_mid) / 2, 0.5, f"{sign}{diff:.1f} pp",
            ha="center", va="center", fontsize=9, fontweight="bold", color="#37474f",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#37474f", lw=0.8, alpha=0.9), zorder=4)

    ax.set_yticks([0, 1])
    ax.set_yticklabels(["After", "Before"])
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of full BUSCO set")
    
    handles = [Patch(facecolor=BUSCO_PAL[k], label=lab) for k, lab in seg]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.45),
              ncol=4, frameon=False, fontsize=9)
# Main figure
def fig_main(master, busco_res, outdir):
    fig = plt.figure(figsize=(15, 9.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)

    ax = fig.add_subplot(gs[0, 0])
    scatter_by_call(ax, master, "median_ks", "depth_ratio_min", KS_SPLIT, ARTI_DEPTH)
    ax.set_xlabel("Median Ks"); ax.set_ylabel("Min depth / median")
    ax.set_title("A.", loc="left", fontweight="bold"); ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[0, 1])
    scatter_by_call(ax, master, "median_ks", "samepos_frac", KS_SPLIT, ARTI_SAMEPOS)
    ax.set_xlabel("Median Ks"); ax.set_ylabel(SHARED_LABEL)
    ax.set_title("B.", loc="left", fontweight="bold")

    draw_decision_tree(fig.add_subplot(gs[0, 2]), leaf_counts(master))
    draw_upset_split(fig, gs[1, 0], master)
    draw_busco_bar(fig.add_subplot(gs[1, 1:3]), busco_res)

    fig.savefig(outdir / "fig_main.pdf", bbox_inches="tight")
    plt.close(fig)

# Supplementary
def fig_supplementary(master, outdir):
    fig, ax = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)

    a = ax[0, 0]
    for call, col in CALL_PAL.items():
        s = master[master["final_call"] == call]
        a.scatter(s["depth_ratio_q"], s["depth_ratio_t"], s=s["n_pairs"] * 3, c=col,
                  alpha=0.7, edgecolor="white", lw=0.5, label=call)
    a.plot([0, 1.5], [0, 1.5], color="gray", ls="--", lw=1)
    a.set_xlabel("Depth q"); a.set_ylabel("Depth t")
    a.set_title("A.", loc="left", fontweight="bold"); a.legend(frameon=False)

    a = ax[0, 1]
    sns.boxplot(data=master, x="label_ks", y="snp_per_kb", hue="loc", hue_order=["intra", "inter"],
                palette=LOC_PAL, ax=a, fliersize=3)
    a.set_title("B.", loc="left", fontweight="bold"); a.set_xlabel("")

    a = ax[1, 0]
    sub = master[master.get("samepos_either", 0) >= MIN_EITHER] if "samepos_either" in master else master
    sns.boxplot(data=sub, x="label_ks", y="samepos_frac", hue="loc", hue_order=["intra", "inter"],
                palette=LOC_PAL, ax=a, fliersize=3)
    a.set_title("C.", loc="left", fontweight="bold"); a.set_xlabel(""); a.set_ylabel(SHARED_LABEL)

    a = ax[1, 1]
    if "category" in master:
        ct = pd.crosstab(master["final_call"], master["category"])
        ct.plot(kind="bar", stacked=True, ax=a,
                color={"private": COL_DUP, "lolium": COL_FRAG, "ancient": COL_COMPLETE})
        a.set_title("D.", loc="left", fontweight="bold")
        a.set_xlabel(""); a.tick_params(axis="x", rotation=0)
    else:
        a.axis("off")

    fig.savefig(outdir / "fig_supplementary.pdf", bbox_inches="tight")
    plt.close(fig)

# Genome landscape
def genome_layout(fai):
    ch = pd.read_csv(fai, sep="\t", header=None).iloc[:, :2]
    ch.columns = ["chr", "len"]
    return ch[ch["chr"].isin(CHROMS)].set_index("chr").reindex(CHROMS).dropna()

def _filled(ax, x, y, base, ylim, clen):
    ax.fill_between(x, base, y, where=(y >= base), facecolor=COL_HIGH, alpha=0.4, interpolate=True)
    ax.fill_between(x, base, y, where=(y < base), facecolor=COL_LOW, alpha=0.4, interpolate=True)
    ax.plot(x, y, color="#333333", lw=0.5)
    ax.plot([0, clen], [base, base], color="gray", lw=0.6, ls="--")  # stops at chr boundary
    ax.set_ylim(*ylim)

def fig_genome_landscape(master, winfis, windepth, masked, outdir, busco_loci=None):
    ch = genome_layout(CONFIG["fai"]); max_mb = ch["len"].max() / 1e6
    bl = busco_loci if busco_loci is not None else load_busco_loci()
    
    fig = plt.figure(figsize=(16, 22.6))
    outer = fig.add_gridspec(len(ch), 1, hspace=0.4)

    for i, chrom in enumerate(ch.index):
        clen = ch.loc[chrom, "len"] / 1e6
        
        # Nested grids shenanigans
        inner_main = outer[i].subgridspec(2, 1, height_ratios=[14, 5], hspace=0.3)
        inner_top = inner_main[0].subgridspec(2, 1, height_ratios=[2, 2], hspace=0.5)
        inner_bot = inner_main[1].subgridspec(3, 1, height_ratios=[2, 2, 2], hspace=0.1)
        A = fig.add_subplot(inner_top[0])
        B = fig.add_subplot(inner_top[1])
        C = fig.add_subplot(inner_bot[0])
        D = fig.add_subplot(inner_bot[1])
        E = fig.add_subplot(inner_bot[2])

        if winfis is not None:
            w = winfis[winfis["chr"] == chrom].sort_values("win")
            _filled(A, w["win"] / 1e6, w["fis"].rolling(SMOOTH, center=True, min_periods=1).mean(), 0, (-1, 1), clen)
        if windepth is not None:
            d = windepth[windepth["chr"] == chrom].sort_values("win")
            _filled(B, d["win"] / 1e6, d["depth_ratio"].rolling(SMOOTH, center=True, min_periods=1).mean(), 1, (0, 3), clen)
            B.plot([0, clen], [COLLAPSE_DEPTH, COLLAPSE_DEPTH], color=COL_HIGH, ls=":", lw=1)

        for tk in (C, D, E):
            tk.add_patch(Rectangle((0, 0.1), clen, 0.8, facecolor="#f0f0f0", edgecolor="#cccccc", lw=0.4, zorder=0))
            tk.set_ylim(0, 1); tk.set_yticks([])

        if bl is not None:
            for kind, col in [("single", COL_COMPLETE), ("dup", COL_DUP)]:
                s = bl[(bl["kind"] == kind) & (bl["chr"] == chrom)]
                C.vlines((s["start"] + s["end"]) / 2 / 1e6, 0.1, 0.9, color=col, lw=0.5, alpha=0.7)

        for _, b in master.iterrows():
            col = CALL_PAL.get(b["final_call"], "gray")
            for c, st, en in [(b["q_chr"], b["q_start"], b["q_end"]), (b["t_chr"], b["t_start"], b["t_end"])]:
                if c == chrom:
                    D.broken_barh([(st / 1e6, (en - st) / 1e6)], (0.1, 0.8), facecolors=col, alpha=0.9)

        for _, m in masked[masked["chr"] == chrom].iterrows():
            E.broken_barh([(m["start"] / 1e6, (m["end"] - m["start"]) / 1e6)], (0.1, 0.8), facecolors=COL_MASK, alpha=0.9)

        for ax, name in zip([A, B, C, D, E], ["FIS", "Depth", "BUSCO", "Call", "Mask"]):
            ax.set_xlim(0, max_mb)
            ax.set_ylabel(name, rotation=0, ha="right", va="center", fontsize=8)
            for sp in ("top", "right", "left"):
                ax.spines[sp].set_visible(False)
            if ax is not E:
                ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
        A.set_title(f"{chrom}", loc="left", fontweight="bold", fontsize=11)
        E.set_xlabel("Position (Mb)", fontsize=8)

    fig.savefig(outdir / "fig_genome_landscape.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--skip-coverage", action="store_true")
    a = ap.parse_args()
    outdir = Path(CONFIG["outdir"]); outdir.mkdir(parents=True, exist_ok=True)

    log("[1/8] Blocks")
    genes = load_genes(); blocks = load_blocks(genes); log(f"  {len(blocks)} blocks")
    master = blocks.copy(); windepth = None

    if not a.skip_coverage:
        log("[2/8] Coverage")
        cov, windepth = coverage_metrics(blocks, outdir, a.force)
        master = master.merge(cov, on="block_id", how="left")
    else:
        log("[2/8] Coverage skipped")

    log("[3/8] Genic SNPs")
    master = master.merge(cached(outdir / "blocks_snp.tsv", lambda: snp_in_genes(blocks, genes), a.force),
                          on="block_id", how="left")

    log("[4/8] FIS")
    hwe = load_hwe()
    master = master.merge(cached(outdir / "blocks_fis.tsv", lambda: fis_per_block(blocks, hwe), a.force),
                          on="block_id", how="left")
    winfis = cached(outdir / "fis_win.tsv", lambda: fis_windows(hwe), a.force)

    log("[5/8] Shared polymorphism")
    try:
        master = master.merge(cached(outdir / "blocks_samepos.tsv", lambda: same_position_metric(blocks), a.force),
                              on="block_id", how="left")
    except Exception as exc:
        log(f"  shared_polymorphism ignored ({exc})")

    log("[6/8] Final call + dupshare")
    master["final_call"] = master.apply(_final_call, axis=1)
    sharing = load_sharing()
    if sharing is not None:
        master = master.merge(sharing, on="block_id", how="left")
    master.to_csv(outdir / "block_master.tsv", sep="\t", index=False)
    log("  " + master["final_call"].value_counts().to_string().replace("\n", "\n  "))

    log("[7/8] Masking + BUSCO")
    busco_loci = load_busco_loci()
    masked = masked_intervals(master, busco_loci)
    
    bed = masked.copy(); bed["name"] = "artifact_block_" + bed["block_id"].astype(str)
    bed[["chr", "start", "end", "name"]].to_csv(outdir / "masked_intervals.bed", sep="\t", index=False, header=False)
    log(f"  {len(masked)} copies to mask -> masked_intervals.bed")
    
    busco_res = busco_recompute(masked)
    if busco_res:
        log(f"  BUSCO Complete {busco_res['S_before']:.1f}% -> {busco_res['S_after']:.1f}% | "
            f"Dup {busco_res['D_before']:.1f}% -> {busco_res['D_after']:.1f}%")

    log("[8/8] Figures")
    if not a.skip_coverage:
        fig_main(master, busco_res, outdir)
        fig_supplementary(master, outdir)
    fig_genome_landscape(master, winfis, windepth, masked, outdir, busco_loci)
    log(f"Finished -> {outdir}/")

if __name__ == "__main__":
    main()