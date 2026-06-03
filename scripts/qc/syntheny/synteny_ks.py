"""Per-block median Ks (NG86) with duplicated-BUSCO overlay. Tremona only.

1. Load CDS, translate to proteins.
2. Parse MCScanX collinearity for anchor pairs.
3. Per pair: global protein alignment -> codon back-translation.
4. One AXT file -> KaKs_Calculator (NG method).
5. Per-block median Ks, converted to age with mu_s.
6. Overlay duplicated BUSCOs onto block coordinates per Ks bin.
7. Plot intra vs inter histogram + cumulative BUSCO curve.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from scipy.interpolate import PchipInterpolator

MU_S      = 5.76174e-9     # De La Torre et al. 2017
HET_KS    = 0.0327         # genome heterozygosity (~3.27%)
MAX_KS    = 1.15
BIN_WIDTH = 0.018
KAKS_BIN  = "KaKs_Calculator"

INTER_COLOR     = "#1976d2"
INTRA_COLOR     = "#9c27b0"
BUSCO_CUM_COLOR = "#ff7f0e"

CODON_TABLE = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F',
    'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y',
    'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C',
    'TGA': '*',
    'TGG': 'W',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
    'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N',
    'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S',
    'AGA': 'R', 'AGG': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D',
    'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
}

HEADER_RE = re.compile(r"## Alignment\s+(\d+):.*?(\S+)&(\S+)\s+(plus|minus)")


def translate(cds):
    aa = []
    for i in range(0, len(cds) - 2, 3):
        x = CODON_TABLE.get(cds[i:i+3])
        if x is None or x == '*':
            break
        aa.append(x)
    return ''.join(aa)


def backtranslate(ap1, ap2, c1, c2):
    out1, out2 = [], []
    i1 = i2 = 0
    for a1, a2 in zip(ap1, ap2):
        out1.append('---' if a1 == '-' else c1[i1*3:i1*3+3])
        if a1 != '-':
            i1 += 1
        out2.append('---' if a2 == '-' else c2[i2*3:i2*3+3])
        if a2 != '-':
            i2 += 1
    return ''.join(out1), ''.join(out2)


def parse_collinearity(path):
    blocks, cur = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("## Alignment"):
                m = HEADER_RE.search(line)
                if m:
                    if cur:
                        blocks.append(cur)
                    cur = {"id": int(m.group(1)), "pairs": []}
            elif cur and line and not line.startswith("#"):
                p = line.split("\t")
                if len(p) >= 3:
                    cur["pairs"].append((p[1].strip(), p[2].strip()))
    if cur:
        blocks.append(cur)
    return blocks


def load_gff(gff_path):
    """gene/mRNA -> (chrom, start, end), keyed by canonical ID."""
    genes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip().split("\t")
            if len(p) == 9 and p[2] in ("gene", "mRNA"):
                m = re.search(r"ID=([^;]+)", p[8])
                if m:
                    gid = re.sub(r"^(rna-|gene-)", "", m.group(1))
                    genes[gid] = (p[0], int(p[3]), int(p[4]))
    return genes


def build_axt(blocks, proteins, cds_seqs, gene_coords, aligner, axt_path):
    """Write all valid pairs to one AXT. Track block extent and real
    chromosome name from the GFF (so BUSCO coordinates line up)."""
    pair_info = {}
    block_coords = {}
    with open(axt_path, "w") as out:
        for b in blocks:
            bc = block_coords.setdefault(b["id"], {
                "q_chr": "?", "t_chr": "?",
                "q_starts": [], "q_ends": [], "t_starts": [], "t_ends": [],
            })
            for g1, g2 in b["pairs"]:
                if g1 in gene_coords:
                    bc["q_chr"] = gene_coords[g1][0]
                    bc["q_starts"].append(gene_coords[g1][1])
                    bc["q_ends"].append(gene_coords[g1][2])
                if g2 in gene_coords:
                    bc["t_chr"] = gene_coords[g2][0]
                    bc["t_starts"].append(gene_coords[g2][1])
                    bc["t_ends"].append(gene_coords[g2][2])

                if g1 not in proteins or g2 not in proteins:
                    continue
                try:
                    aln = aligner.align(proteins[g1], proteins[g2])[0]
                    ca1, ca2 = backtranslate(str(aln[0]), str(aln[1]),
                                             cds_seqs[g1], cds_seqs[g2])
                except Exception:
                    continue
                if len(ca1) < 90:
                    continue
                name = f"B{b['id']}__{g1}__{g2}"
                pair_info[name] = b["id"]
                out.write(f"{name}\n{ca1}\n{ca2}\n\n")
    return pair_info, block_coords


def run_kakscalculator(axt_path, out_path):
    cmd = [KAKS_BIN, "-i", axt_path, "-o", out_path, "-c", "1", "-m", "NG"]
    print(f"[info] running: {' '.join(cmd)}", file=sys.stderr)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.stderr.strip():
        print("[kaks stderr]\n" + res.stderr, file=sys.stderr)
    if res.returncode != 0:
        raise RuntimeError(f"KaKs_Calculator exited with code {res.returncode}")
    if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
        raise RuntimeError(
            f"KaKs_Calculator wrote no output to {out_path}. "
            "Rerun with --keep-tmp and inspect the AXT file."
        )


def parse_kaks_output(path, pair_info):
    DEBUG_BLOCK = 46          # mettre None pour couper les messages
    results = []
    with open(path) as fh:
        header = fh.readline().rstrip().split("\t")
        i_seq = header.index("Sequence")
        i_ka  = header.index("Ka")
        i_ks  = header.index("Ks")
        i_sub = header.index("Substitutions")
        if DEBUG_BLOCK is not None:
            print(f"[dbg] indices Ka={i_ka} Ks={i_ks} Substitutions={i_sub}",
                  file=sys.stderr)
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) <= max(i_seq, i_ka, i_ks, i_sub):
                continue
            name = p[i_seq]
            if name not in pair_info:
                continue
            bid = pair_info[name]
            dbg = (DEBUG_BLOCK is not None and bid == DEBUG_BLOCK)
            try:
                ka, ks = float(p[i_ka]), float(p[i_ks])
                branch = "float"
            except ValueError:
                try:
                    sub = float(p[i_sub])
                except ValueError:
                    sub = None
                if sub == 0:
                    ka, ks, branch = 0.0, 0.0, "identical->0"
                else:
                    if dbg:
                        print(f"[dbg] {name}: Ks={p[i_ks]!r} "
                              f"sub={p[i_sub]!r} -> JETE (NA, sub!=0)",
                              file=sys.stderr)
                    continue
            if not (0 <= ks < 5):
                if dbg:
                    print(f"[dbg] {name}: ks={ks} -> JETE (hors borne)",
                          file=sys.stderr)
                continue
            if dbg:
                print(f"[dbg] {name}: ks={ks} ({branch}) -> GARDE",
                      file=sys.stderr)
            results.append({"block_id": bid, "ka": ka, "ks": ks})
    return results

def aggregate_per_block(pair_results, block_coords):
    by_block = defaultdict(list)
    for r in pair_results:
        by_block[r["block_id"]].append(r)
    out = []
    for bid, pairs in by_block.items():
        bc = block_coords.get(bid, {})
        out.append({
            "id": bid,
            "q_chr": bc.get("q_chr", "?"),
            "t_chr": bc.get("t_chr", "?"),
            "q_start": min(bc["q_starts"]) if bc.get("q_starts") else -1,
            "q_end":   max(bc["q_ends"])   if bc.get("q_ends")   else -1,
            "t_start": min(bc["t_starts"]) if bc.get("t_starts") else -1,
            "t_end":   max(bc["t_ends"])   if bc.get("t_ends")   else -1,
            "n_pairs": len(pairs),
            "median_ks": float(np.median([p["ks"] for p in pairs])),
            "median_ka": float(np.median([p["ka"] for p in pairs])),
        })
    return out


def main(args):
    print("[info] loading CDS", file=sys.stderr)
    cds_seqs, proteins = {}, {}
    for rec in SeqIO.parse(args.cds, "fasta"):
        c = str(rec.seq).upper()
        c = c[:len(c) - len(c) % 3]
        p = translate(c)
        if len(p) >= 30:
            cds_seqs[rec.id] = c[:len(p) * 3]
            proteins[rec.id] = p

    blocks = parse_collinearity(args.collinearity)
    print(f"[info] {len(blocks)} synteny blocks", file=sys.stderr)

    gene_coords = load_gff(args.gff)

    busco_df = pd.read_csv(
        args.busco, sep='\t', comment='#',
        names=['Busco_id', 'Status', 'Sequence',
               'Gene_Start', 'Gene_End', 'Strand', 'Score', 'Length'],
    )
    busco_df['Gene_Start'] = pd.to_numeric(busco_df['Gene_Start'], errors='coerce')
    busco_df['Gene_End']   = pd.to_numeric(busco_df['Gene_End'],   errors='coerce')
    dup_buscos = busco_df[busco_df['Status'] == 'Duplicated'].dropna(
        subset=['Gene_Start', 'Gene_End']
    )
    total_dup = dup_buscos['Busco_id'].nunique()

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    out_dir = os.path.dirname(os.path.abspath(args.output)) or "."
    os.makedirs(out_dir, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="kaks_", dir=out_dir)
    print(f"[info] work dir: {tmp}", file=sys.stderr)
    try:
        axt = os.path.join(tmp, "pairs.axt")
        kaks = os.path.join(tmp, "pairs.kaks")
        print("[info] aligning pairs and building AXT", file=sys.stderr)
        pair_info, block_coords = build_axt(
            blocks, proteins, cds_seqs, gene_coords, aligner, axt
        )
        print(f"[info] {len(pair_info)} pairs prepared", file=sys.stderr)
        run_kakscalculator(axt, kaks)
        pair_results = parse_kaks_output(kaks, pair_info)
        print(f"[info] {len(pair_results)} pairs with valid Ks", file=sys.stderr)
    finally:
        if args.keep_tmp:
            print(f"[info] kept work dir: {tmp}", file=sys.stderr)
        else:
            shutil.rmtree(tmp, ignore_errors=True)

    block_results = aggregate_per_block(pair_results, block_coords)
    n_located = sum(1 for r in block_results if r["q_start"] != -1)
    print(f"[info] {len(block_results)} blocks with Ks, "
          f"{n_located} with GFF coordinates", file=sys.stderr)

    # TSV
    tsv = args.output.replace('.pdf', '.tsv')
    with open(tsv, "w") as fh:
        fh.write("block_id\tq_chr\tt_chr\tq_start\tq_end\tt_start\tt_end"
                 "\tn_pairs\tmedian_ks\tmedian_ka\tage_Mya\n")
        for r in block_results:
            age = r["median_ks"] / (2 * MU_S) / 1e6
            fh.write(f"{r['id']}\t{r['q_chr']}\t{r['t_chr']}"
                     f"\t{r['q_start']}\t{r['q_end']}\t{r['t_start']}\t{r['t_end']}"
                     f"\t{r['n_pairs']}\t{r['median_ks']:.4f}\t{r['median_ka']:.4f}"
                     f"\t{age:.1f}\n")
    print(f"[info] wrote {tsv}", file=sys.stderr)

    intra_ks = [r["median_ks"] for r in block_results if r["q_chr"] == r["t_chr"]]
    inter_ks = [r["median_ks"] for r in block_results if r["q_chr"] != r["t_chr"]]

    bins = np.arange(0, MAX_KS + BIN_WIDTH, BIN_WIDTH)
    centres = (bins[:-1] + bins[1:]) / 2

    seen_cum = set()
    cumulative = []
    for k in range(len(bins) - 1):
        lo, hi = bins[k], bins[k+1]
        for r in block_results:
            if not (lo <= r['median_ks'] < hi):
                continue
            for chrom, s, e in ((r['q_chr'], r['q_start'], r['q_end']),
                                (r['t_chr'], r['t_start'], r['t_end'])):
                if s == -1:
                    continue
                over = dup_buscos[
                    (dup_buscos['Sequence']   == chrom) &
                    (dup_buscos['Gene_Start'] <= e) &
                    (dup_buscos['Gene_End']   >= s)
                ]
                seen_cum.update(over['Busco_id'].tolist())
        cumulative.append(len(seen_cum))
    cumulative = np.array(cumulative, dtype=float)
    cumulative_pct = (cumulative / total_dup * 100) if total_dup > 0 \
        else np.zeros_like(cumulative)

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.labelsize": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "legend.fontsize": 9,
    })
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8.5, 8.5), sharex=True,
        gridspec_kw={'height_ratios': [3, 1.2]},
    )
    plt.subplots_adjust(hspace=0.1)

    inter_c, _ = np.histogram(inter_ks, bins=bins)
    intra_c, _ = np.histogram(intra_ks, bins=bins)
    bw = BIN_WIDTH * 0.95

    ax1.bar(centres, inter_c, width=bw, color=INTER_COLOR, edgecolor="white",
            linewidth=0.2, alpha=0.75,
            label=f"Inter-chromosomal (n = {len(inter_ks)})", zorder=2)
    ax1.bar(centres, intra_c, width=bw, bottom=inter_c, color=INTRA_COLOR,
            edgecolor="white", linewidth=0.2, alpha=0.75,
            label=f"Intra-chromosomal (n = {len(intra_ks)})", zorder=2)

    y_top = max(int((inter_c + intra_c).max()), 1) * 1.05
    ax1.axvline(HET_KS, color="#888", ls="--", lw=1, zorder=10)
    ax1.axvline(5.4 * 2 * MU_S * 1e6, color="#c73838", ls="--", lw=1, zorder=10)
    secax = ax1.secondary_xaxis(
        "top",
        functions=(lambda ks: ks / (2 * MU_S) / 1e6,
                   lambda my: my * 2 * MU_S * 1e6),
    )
    secax.set_xlabel("Estimated divergence time (Mya)")
    ax1.set_ylabel("Number of synteny blocks")
    ax1.set_xlim(0, MAX_KS)
    ax1.set_ylim(0, y_top)
    ax1.legend(loc="upper right")

    x_smooth = np.linspace(centres.min(), centres.max(), 300)
    y_smooth = PchipInterpolator(centres, cumulative_pct)(x_smooth)
    ax2.plot(x_smooth, y_smooth, color=BUSCO_CUM_COLOR, linewidth=2.5)
    ax2.fill_between(x_smooth, 0, y_smooth, color=BUSCO_CUM_COLOR, alpha=0.15)
    ax2.axvline(HET_KS, color="#888", ls="--", lw=1, zorder=10)
    ax2.axvline(5.4 * 2 * MU_S * 1e6, color="#c73838", ls="--", lw=1, zorder=10)
    ax2.set_xlabel("Median Ks per synteny block (NG86)")
    ax2.set_ylabel("% of duplicated\nBUSCOs covered",
                   color=BUSCO_CUM_COLOR, fontsize=10, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor=BUSCO_CUM_COLOR)
    ax2.set_ylim(0, 100)

    plt.savefig(args.output, bbox_inches="tight")
    plt.savefig(args.output.rsplit(".", 1)[0] + ".png", bbox_inches="tight", dpi=300)
    covered = cumulative_pct[-1] if len(cumulative_pct) else 0.0
    print(f"[info] wrote {args.output} "
          f"({covered:.1f}% of {total_dup} duplicated BUSCOs covered by blocks)",
          file=sys.stderr)


if __name__ == "__main__":
    D = "results/synteny/tremona"
    p = argparse.ArgumentParser(description="Per-block Ks (NG) + BUSCO overlay, Tremona.")
    p.add_argument("--collinearity", default=f"{D}/tremona.collinearity")
    p.add_argument("--cds",          default=f"{D}/cds.fa")
    p.add_argument("--gff",          default=f"{D}/annotation.gff3")
    p.add_argument("--busco",        default="reference_data/lmultiflorum.tremona_full_table_busco_format.tsv")
    p.add_argument("--output",       default=f"{D}/synteny_ks.pdf")
    p.add_argument("--keep-tmp",     action="store_true")
    main(p.parse_args())