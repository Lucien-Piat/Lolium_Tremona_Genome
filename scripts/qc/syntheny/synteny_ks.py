"""Per-block median Ks via KaKs_Calculator (YN00), with BUSCO overlay.

Workflow:
1. Load CDS, translate to proteins ourselves (keeps prot/CDS aligned).
2. Parse MCScanX collinearity for anchor pairs.
3. For each pair: protein alignment (BLOSUM62) -> codon back-translation.
4. Batch all pairs into one AXT file, run KaKs_Calculator with YN00 + NG86.
5. Aggregate per-block median Ks. Convert to age using grass mu_s.
6. Map duplicated BUSCOs to physical coordinates of blocks per Ks bin.
7. Plot intra vs inter distribution with cumulative BUSCO overlay.
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

MU_S = 5.76174e-9     # De La Torre et al. 2017
HET_KS = 0.0327       # genome heterozygosity (~3.27%).
MAX_KS = 3.2
BIN_WIDTH = 0.04

KAKS_BIN       = "KaKs_Calculator"
KAKS_METHODS   = "YN,NG"
PRIMARY_METHOD = "YN"

INTER_COLOR     = "#1976d2"
INTRA_COLOR     = "#9c27b0"
BUSCO_CUM_COLOR = "#ff7f0e"

CODON_TABLE = {
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',    # Serine
    'TTC': 'F', 'TTT': 'F',    # Phenylalanine
    'TTA': 'L', 'TTG': 'L',    # Leucine
    'TAC': 'Y', 'TAT': 'Y',    # Tirosine
    'TAA': '*', 'TAG': '*',    # Stop
    'TGC': 'C', 'TGT': 'C',    # Cisteine
    'TGA': '*',    # Stop
    'TGG': 'W',    # Tryptofan
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',    # Leucine
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',    # Proline
    'CAC': 'H', 'CAT': 'H',    # Histidine
    'CAA': 'Q', 'CAG': 'Q',    # Glutamine
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',    # Arginine
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',    # Isoleucine
    'ATG': 'M',    # Methionine
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',    # Threonine
    'AAC': 'N', 'AAT': 'N',    # Asparagine
    'AAA': 'K', 'AAG': 'K',    # Lysine
    'AGC': 'S', 'AGT': 'S',    # Serine
    'AGA': 'R', 'AGG': 'R',    # Arginine
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',    # Valine
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',    # Alanine
    'GAC': 'D', 'GAT': 'D',    # Aspartic Acid
    'GAA': 'E', 'GAG': 'E',    # Glutamic Acid
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'     # Glycine
}

HEADER_RE = re.compile(r"## Alignment\s+(\d+):.*?(\S+)&(\S+)\s+(plus|minus)")


def translate(cds):
    aa = []
    for i in range(0, len(cds) - 2, 3):
        c = cds[i:i+3]
        x = CODON_TABLE.get(c)
        if x is None or x == '*': break
        aa.append(x)
    return ''.join(aa)

def backtranslate(ap1, ap2, c1, c2):
    out1, out2 = [], []
    i1 = i2 = 0
    for a1, a2 in zip(ap1, ap2):
        out1.append('---' if a1 == '-' else c1[i1*3:i1*3+3])
        if a1 != '-': i1 += 1
        out2.append('---' if a2 == '-' else c2[i2*3:i2*3+3])
        if a2 != '-': i2 += 1
    return ''.join(out1), ''.join(out2)

def parse_collinearity(path):
    blocks, cur = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("## Alignment"):
                m = HEADER_RE.search(line)
                if m:
                    if cur: blocks.append(cur)
                    cur = {"id": int(m.group(1)),
                           "q_chr": m.group(2), "t_chr": m.group(3),
                           "pairs": []}
            elif cur and line and not line.startswith("#"):
                p = line.split("\t")
                if len(p) >= 3:
                    cur["pairs"].append((p[1].strip(), p[2].strip()))
    if cur: blocks.append(cur)
    return blocks


def load_gff(gff_path):
    genes = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.rstrip().split("\t")
            if len(p) == 9 and p[2] in ["gene", "mRNA"]:
                chrom, start, end = p[0], int(p[3]), int(p[4])
                match = re.search(r"ID=([^;]+)", p[8])
                if match: genes[match.group(1)] = (chrom, start, end)
    return genes

def age_to_ks(age_mya):
    return age_mya * 2 * MU_S * 1e6

def build_axt(blocks, proteins, cds_seqs, gene_coords, aligner, axt_path):
    """
    Write a single AXT file with all valid pairs.
    Returns (pair_info, block_coords).
    Block coordinates are tracked from ALL anchor genes, regardless of
    whether the alignment for Ks computation succeeded.
    """
    pair_info = {}
    block_coords = {}

    with open(axt_path, "w") as out:
        for b in blocks:
            bc = block_coords.setdefault(b["id"], {
                "q_chr": b["q_chr"], "t_chr": b["t_chr"],
                "q_starts": [], "q_ends": [],
                "t_starts": [], "t_ends": [],
            })
            for g1, g2 in b["pairs"]:
                # Track block extent from gene coords regardless of alignment
                if g1 in gene_coords:
                    bc["q_starts"].append(gene_coords[g1][1])
                    bc["q_ends"].append(gene_coords[g1][2])
                if g2 in gene_coords:
                    bc["t_starts"].append(gene_coords[g2][1])
                    bc["t_ends"].append(gene_coords[g2][2])

                if g1 not in proteins or g2 not in proteins:
                    continue
                try:
                    aln = aligner.align(proteins[g1], proteins[g2])[0]
                    ap1, ap2 = str(aln[0]), str(aln[1])
                    ca1, ca2 = backtranslate(ap1, ap2, cds_seqs[g1], cds_seqs[g2])
                except Exception:
                    continue
                if len(ca1) < 90:
                    continue
                name = f"B{b['id']}__{g1}__{g2}"
                pair_info[name] = (b["id"], b["q_chr"], b["t_chr"], g1, g2)
                out.write(f"{name}\n{ca1}\n{ca2}\n\n")
    return pair_info, block_coords


def run_kakscalculator(axt_path, out_path):
    cmd = [KAKS_BIN, "-i", axt_path, "-o", out_path, "-c", "1"]
    for method in KAKS_METHODS.split(","):
        method = method.strip()
        if method:
            cmd += ["-m", method]

    print(f"[info] running: {' '.join(cmd)}", file=sys.stderr)
    res = subprocess.run(cmd, capture_output=True, text=True)

    if res.stdout.strip():
        print("[kaks stdout]\n" + res.stdout, file=sys.stderr)
    if res.stderr.strip():
        print("[kaks stderr]\n" + res.stderr, file=sys.stderr)

    if res.returncode != 0:
        raise RuntimeError(f"KaKs_Calculator exited with code {res.returncode}")

    if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
        raise RuntimeError(
            f"KaKs_Calculator returned 0 but wrote no output to {out_path}. "
            "Check the [kaks ...] messages above (often an unrecognised -m "
            "method, an unsupported option, or a malformed AXT record). "
            "Rerun with --keep-tmp to inspect the AXT file."
        )


def parse_kaks_output(path, pair_info):
    """Parse KaKs_Calculator tabular output.
    Returns {method: [{block_id, q_chr, t_chr, g1, g2, ka, ks}, ...]}.
    """
    results = defaultdict(list)
    with open(path) as fh:
        header = fh.readline().rstrip().split("\t")
        try:
            i_seq    = header.index("Sequence")
            i_method = header.index("Method")
            i_ka     = header.index("Ka")
            i_ks     = header.index("Ks")
        except ValueError:
            print(f"[error] unexpected KaKs header: {header}", file=sys.stderr)
            raise
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) <= max(i_seq, i_method, i_ka, i_ks):
                continue
            name = p[i_seq]; method = p[i_method]
            try:
                ka = float(p[i_ka]); ks = float(p[i_ks])
            except ValueError:
                continue
            if not (0 < ks < 5):
                continue
            if name not in pair_info:
                continue
            block_id, q_chr, t_chr, g1, g2 = pair_info[name]
            results[method].append({
                "block_id": block_id, "q_chr": q_chr, "t_chr": t_chr,
                "g1": g1, "g2": g2, "ka": ka, "ks": ks,
            })
    return results


def aggregate_per_block(pair_results, block_coords):
    by_block = defaultdict(list)
    for r in pair_results:
        by_block[r["block_id"]].append(r)
    out = []
    for bid, pairs in by_block.items():
        bc = block_coords.get(bid, {})
        q_s = min(bc["q_starts"]) if bc.get("q_starts") else -1
        q_e = max(bc["q_ends"])   if bc.get("q_ends")   else -1
        t_s = min(bc["t_starts"]) if bc.get("t_starts") else -1
        t_e = max(bc["t_ends"])   if bc.get("t_ends")   else -1
        ks_vals = [p["ks"] for p in pairs]
        ka_vals = [p["ka"] for p in pairs]
        out.append({
            "id": bid,
            "q_chr": pairs[0]["q_chr"], "t_chr": pairs[0]["t_chr"],
            "q_start": q_s, "q_end": q_e,
            "t_start": t_s, "t_end": t_e,
            "n_pairs": len(pairs),
            "median_ks": float(np.median(ks_vals)),
            "median_ka": float(np.median(ka_vals)),
        })
    return out


def main(args):
    print("[info] loading CDS", file=sys.stderr)
    cds_seqs, proteins = {}, {}
    for rec in SeqIO.parse(args.cds, "fasta"):
        c = str(rec.seq).upper()
        c = c[:len(c) - len(c) % 3]
        p = translate(c)
        c = c[:len(p) * 3]
        if len(p) >= 30:
            cds_seqs[rec.id] = c
            proteins[rec.id] = p

    blocks = parse_collinearity(args.collinearity)
    print(f"[info] {len(blocks)} synteny blocks", file=sys.stderr)

    print(f"[info] loading GFF and BUSCO data", file=sys.stderr)
    gene_coords = load_gff(args.gff)

    busco_df = pd.read_csv(
        args.busco, sep='\t', comment='#',
        names=['Busco_id', 'Status', 'Sequence',
               'Gene_Start', 'Gene_End', 'Strand', 'Score', 'Length'],
    )
    busco_df['Sequence']   = busco_df['Sequence'].str.replace('chr', 'lm')
    busco_df['Gene_Start'] = pd.to_numeric(busco_df['Gene_Start'], errors='coerce')
    busco_df['Gene_End']   = pd.to_numeric(busco_df['Gene_End'],   errors='coerce')
    dup_buscos = busco_df[busco_df['Status'] == 'Duplicated'].dropna(
        subset=['Gene_Start', 'Gene_End']
    )
    total_unique_dup_buscos = dup_buscos['Busco_id'].nunique()
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    # Create the temp dir INSIDE the Singularity bind 
    out_dir = os.path.dirname(os.path.abspath(args.output)) or "."
    base_tmp = args.workdir or out_dir
    os.makedirs(base_tmp, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="kaks_", dir=base_tmp)
    print(f"[info] work dir: {tmp}", file=sys.stderr)

    try:
        axt = os.path.join(tmp, "pairs.axt")
        out = os.path.join(tmp, "pairs.kaks")

        print("[info] aligning pairs and building AXT", file=sys.stderr)
        pair_info, block_coords = build_axt(
            blocks, proteins, cds_seqs, gene_coords, aligner, axt
        )
        print(f"[info] {len(pair_info)} pairs prepared", file=sys.stderr)

        run_kakscalculator(axt, out)

        results = parse_kaks_output(out, pair_info)
        for m, lst in results.items():
            print(f"[info] {m}: {len(lst)} pairs with valid Ks",
                  file=sys.stderr)
    finally:
        if args.keep_tmp:
            print(f"[info] kept work dir: {tmp}", file=sys.stderr)
        else:
            shutil.rmtree(tmp, ignore_errors=True)

    # Primary aggregation
    block_results = aggregate_per_block(
        results.get(PRIMARY_METHOD, []), block_coords
    )
    print(f"[info] {len(block_results)} blocks with valid Ks", file=sys.stderr)

    ng_blocks = aggregate_per_block(results.get("NG", []), block_coords)
    ng_by_id = {b["id"]: b for b in ng_blocks}

    # Write TSV with both methods
    tsv = args.output.replace('.pdf', '.tsv')
    with open(tsv, "w") as out:
        out.write("block_id\tq_chr\tt_chr\tq_start\tq_end\tt_start\tt_end"
                  "\tn_pairs\tmedian_ks_yn\tmedian_ka_yn"
                  "\tmedian_ks_ng\tmedian_ka_ng\tage_Mya_yn\n")
        for r in block_results:
            ng = ng_by_id.get(r["id"])
            ng_ks = "NA" if ng is None else f"{ng['median_ks']:.4f}"
            ng_ka = "NA" if ng is None else f"{ng['median_ka']:.4f}"
            age   = r["median_ks"] / (2 * MU_S) / 1e6
            out.write(
                f"{r['id']}\t{r['q_chr']}\t{r['t_chr']}"
                f"\t{r['q_start']}\t{r['q_end']}\t{r['t_start']}\t{r['t_end']}"
                f"\t{r['n_pairs']}"
                f"\t{r['median_ks']:.4f}\t{r['median_ka']:.4f}"
                f"\t{ng_ks}\t{ng_ka}\t{age:.1f}\n"
            )
    print(f"[info] wrote {tsv}", file=sys.stderr)

    intra_ks = [r["median_ks"] for r in block_results if r["q_chr"] == r["t_chr"]]
    inter_ks = [r["median_ks"] for r in block_results if r["q_chr"] != r["t_chr"]]

    bins = np.arange(0, MAX_KS + BIN_WIDTH, BIN_WIDTH)
    centres = (bins[:-1] + bins[1:]) / 2
    busco_counts = []

    for k in range(len(bins) - 1):
        bin_min, bin_max = bins[k], bins[k+1]
        blocks_in_bin = [r for r in block_results
                         if bin_min <= r['median_ks'] < bin_max]

        unique_buscos = set()
        for r in blocks_in_bin:
            if r['q_start'] != -1:
                over_q = dup_buscos[
                    (dup_buscos['Sequence']   == r['q_chr']) &
                    (dup_buscos['Gene_Start'] <= r['q_end']) &
                    (dup_buscos['Gene_End']   >= r['q_start'])
                ]
                unique_buscos.update(over_q['Busco_id'].tolist())
            if r['t_start'] != -1:
                over_t = dup_buscos[
                    (dup_buscos['Sequence']   == r['t_chr']) &
                    (dup_buscos['Gene_Start'] <= r['t_end']) &
                    (dup_buscos['Gene_End']   >= r['t_start'])
                ]
                unique_buscos.update(over_t['Busco_id'].tolist())
        busco_counts.append(len(unique_buscos))

    cumulative_buscos = np.cumsum(busco_counts)
    if total_unique_dup_buscos > 0:
        cumulative_pct = (cumulative_buscos / total_unique_dup_buscos) * 100
    else:
        cumulative_pct = np.zeros_like(cumulative_buscos)

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

    # TOP: histogram
    ax1.bar(centres, inter_c, width=bw, color=INTER_COLOR,
            edgecolor="white", linewidth=0.2, alpha=0.75,
            label=f"Inter-chromosomal (n = {len(inter_ks)})", zorder=2)
    ax1.bar(centres, intra_c, width=bw, color=INTRA_COLOR,
            bottom=inter_c, edgecolor="white", linewidth=0.2, alpha=0.75,
            label=f"Intra-chromosomal (n = {len(intra_ks)})", zorder=2)


    y_top = max(int((inter_c + intra_c).max()), 1) * 1.05
    ax1.axvline(HET_KS, color="#888", ls="--", lw=1, zorder=10)
    ax1.axvline(age_to_ks(5.4), color="#c73838", ls="--", lw=1, zorder=10)

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

    # BOTTOM: cumulative BUSCO curve
    x_smooth = np.linspace(centres.min(), centres.max(), 300)
    pchip = PchipInterpolator(centres, cumulative_buscos)
    y_smooth = pchip(x_smooth)

    ax2.plot(x_smooth, y_smooth, color=BUSCO_CUM_COLOR, linewidth=2.5)
    ax2.fill_between(x_smooth, 0, y_smooth, color=BUSCO_CUM_COLOR, alpha=0.15)

    ax2.axvline(HET_KS, color="#888", ls="--", lw=1, zorder=10)
    ax2.axvline(age_to_ks(5.4), color="#c73838", ls="--", lw=1, zorder=10)
    
    ax2.set_xlabel(f"Median Ks per synteny block (NG86)")
    ax2.set_ylabel("Cumulative\nDuplicated BUSCOs",
                   color=BUSCO_CUM_COLOR, fontsize=10, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor=BUSCO_CUM_COLOR)
    ax2.set_ylim(0, max(cumulative_buscos) * 1.1 if max(cumulative_buscos) > 0 else 1)

    plt.savefig(args.output, bbox_inches="tight")
    png = args.output.rsplit(".", 1)[0] + ".png"
    plt.savefig(png, bbox_inches="tight", dpi=300)
    print(f"[info] wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--collinearity", required=True)
    p.add_argument("--cds",          required=True)
    p.add_argument("--gff",          default="results/data_circo/tremona.gene_annotation.placed.gff")
    p.add_argument("--busco",        default="results/data_circo/full_table_busco_format.tsv")
    p.add_argument("--output",       default="results/data_circo/synteny_ks.pdf")
    p.add_argument("--workdir",      default=None,
                   help="Directory for intermediate AXT/kaks files. Defaults "
                        "to the --output directory. Keep it inside the "
                        "Singularity bind path.")
    p.add_argument("--keep-tmp",     action="store_true",
                   help="Do not delete the work dir (useful for debugging).")
    main(p.parse_args())