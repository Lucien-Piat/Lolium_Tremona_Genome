import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore
from Bio import SeqIO # type: ignore
from Bio.Align import PairwiseAligner, substitution_matrices # type: ignore

MAX_KS    = 3.2
BIN_WIDTH = 0.04
KAKS_BIN  = "KaKs_Calculator"

INTER_COLOR = "#1976d2"
INTRA_COLOR = "#9c27b0"

SYN_BASE = "results/synteny"
DATASETS = ["tremona", "rabiosa","paraquat", "perenne", "brachypodium, oryza"]
OUTPUT   = "results/synteny/synteny_ks_all.pdf"

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
    """Return [{id, q_chr, t_chr, pairs}]. q_chr/t_chr come from the header."""
    blocks, cur = [], None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("## Alignment"):
                m = HEADER_RE.search(line)
                if m:
                    if cur:
                        blocks.append(cur)
                    cur = {"id": int(m.group(1)),
                           "q_chr": m.group(2), "t_chr": m.group(3), "pairs": []}
            elif cur and line and not line.startswith("#"):
                p = line.split("\t")
                if len(p) >= 3:
                    cur["pairs"].append((p[1].strip(), p[2].strip()))
    if cur:
        blocks.append(cur)
    return blocks


def build_axt(blocks, proteins, cds_seqs, aligner, axt_path):
    pair_info = {}
    with open(axt_path, "w") as out:
        for b in blocks:
            for g1, g2 in b["pairs"]:
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
    return pair_info


def run_kakscalculator(axt_path, out_path):
    cmd = [KAKS_BIN, "-i", axt_path, "-o", out_path, "-c", "1", "-m", "NG"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.stderr.strip():
        print("[kaks stderr]\n" + res.stderr, file=sys.stderr)
    if res.returncode != 0:
        raise RuntimeError(f"KaKs_Calculator exited with code {res.returncode}")
    if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
        raise RuntimeError(f"KaKs_Calculator wrote no output to {out_path}.")


def parse_kaks_output(path, pair_info):
    results = []
    with open(path) as fh:
        header = fh.readline().rstrip().split("\t")
        i_seq = header.index("Sequence")
        i_ks  = header.index("Ks")
        i_sub = header.index("Substitutions")
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) <= max(i_seq, i_ks, i_sub) or p[i_seq] not in pair_info:
                continue
            try:
                ks = float(p[i_ks])
            except ValueError:
                try:
                    ks = 0.0 if float(p[i_sub]) == 0 else None
                except ValueError:
                    ks = None
                if ks is None:
                    continue
            if not (0 <= ks < 5):
                continue
            results.append((pair_info[p[i_seq]], ks))
    return results


def compute_genome(name, aligner, tmp_base):
    """Return (inter_ks_list, intra_ks_list) of per-block median Ks."""
    syn = os.path.join(SYN_BASE, name)
    cds_seqs, proteins = {}, {}
    for rec in SeqIO.parse(os.path.join(syn, "cds.fa"), "fasta"):
        c = str(rec.seq).upper()
        c = c[:len(c) - len(c) % 3]
        p = translate(c)
        if len(p) >= 30:
            cds_seqs[rec.id] = c[:len(p) * 3]
            proteins[rec.id] = p

    blocks = parse_collinearity(os.path.join(syn, f"{name}.collinearity"))
    block_chr = {b["id"]: (b["q_chr"], b["t_chr"]) for b in blocks}

    tmp = tempfile.mkdtemp(prefix=f"kaks_{name}_", dir=tmp_base)
    try:
        axt = os.path.join(tmp, "pairs.axt")
        out = os.path.join(tmp, "pairs.kaks")
        pair_info = build_axt(blocks, proteins, cds_seqs, aligner, axt)
        run_kakscalculator(axt, out)
        pairs = parse_kaks_output(out, pair_info)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    by_block = defaultdict(list)
    for bid, ks in pairs:
        by_block[bid].append(ks)

    inter, intra = [], []
    for bid, ks_list in by_block.items():
        med = float(np.median(ks_list))
        qc, tc = block_chr.get(bid, ("?", "?"))
        (intra if qc == tc else inter).append(med)
    return inter, intra


def main():
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    out_dir = os.path.dirname(os.path.abspath(OUTPUT)) or "."
    os.makedirs(out_dir, exist_ok=True)

    data = {}
    for name in DATASETS:
        col = os.path.join(SYN_BASE, name, f"{name}.collinearity")
        cds = os.path.join(SYN_BASE, name, "cds.fa")
        if not (os.path.isfile(col) and os.path.isfile(cds)):
            print(f"[warn] skipping {name}: missing collinearity or cds.fa",
                  file=sys.stderr)
            continue
        print(f"[info] {name}", file=sys.stderr)
        inter, intra = compute_genome(name, aligner, out_dir)
        data[name] = (inter, intra)
        print(f"[info]   {len(inter)} inter, {len(intra)} intra blocks",
              file=sys.stderr)

    if not data:
        sys.exit("[error] no genome produced results")

    bins = np.arange(0, MAX_KS + BIN_WIDTH, BIN_WIDTH)
    centres = (bins[:-1] + bins[1:]) / 2
    bw = BIN_WIDTH * 0.95

    # one shared y limit across all genomes
    hist, ymax = {}, 1
    for name, (inter, intra) in data.items():
        ic, _ = np.histogram(inter, bins=bins)
        ac, _ = np.histogram(intra, bins=bins)
        hist[name] = (ic, ac)
        ymax = max(ymax, int((ic + ac).max()))
    ymax *= 1.05

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.labelsize": 11,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "legend.fontsize": 9,
    })

    n = len(data)
    fig, axes = plt.subplots(n, 1, figsize=(8.5, 1.9 * n + 0.5),
                             sharex=True, sharey=True)
    if n == 1:
        axes = [axes]

    for ax, (name, (inter, intra)) in zip(axes, data.items()):
        ic, ac = hist[name]
        ax.bar(centres, ic, width=bw, color=INTER_COLOR, edgecolor="white",
               linewidth=0.2, alpha=0.8, label=f"Inter (n = {len(inter)})")
        ax.bar(centres, ac, width=bw, bottom=ic, color=INTRA_COLOR,
               edgecolor="white", linewidth=0.2, alpha=0.8,
               label=f"Intra (n = {len(intra)})")
        ax.set_xlim(0, MAX_KS)
        ax.set_ylim(0, ymax)
        ax.text(0.01, 0.90, name, transform=ax.transAxes, ha="left", va="top",
                fontsize=11, fontweight="bold")
        ax.legend(loc="upper right")

    axes[-1].set_xlabel("Median Ks per synteny block (NG86)")
    fig.supylabel("Number of synteny blocks")
    fig.tight_layout()
    fig.savefig(OUTPUT, bbox_inches="tight")
    fig.savefig(OUTPUT.rsplit(".", 1)[0] + ".png", bbox_inches="tight", dpi=300)
    print(f"[info] wrote {OUTPUT}", file=sys.stderr)


if __name__ == "__main__":
    main()