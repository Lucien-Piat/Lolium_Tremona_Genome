import os

import pandas as pd

MIN_FRAC = 0.5
_STR_COLS = ("src_chrom", "dup_chrom", "q_chr", "t_chr")


def read_truth(path):
    # Sequence names must stay str: a chromosome named "1" is otherwise read as
    # int64 and never matches the string keys in a BED.
    df = pd.read_csv(path, sep="\t", dtype={c: str for c in _STR_COLS})
    for c in ("src_start", "src_end", "dup_start", "dup_end", "n_genes"):
        df[c] = df[c].astype(int)
    df["injected_ks"] = df["injected_ks"].astype(float)
    return df


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def read_bed(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) >= 3 and not line.startswith(("#", "track")):
                out.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    return {k: merge(v) for k, v in out.items()}


def covered(bed, chrom, start, end):
    total = 0
    for s, e in bed.get(chrom, ()):
        if e <= start:
            continue
        if s >= end:
            break
        total += min(e, end) - max(s, start)
    return total


def block_detected(bed, row, min_frac=MIN_FRAC):
    # Either copy counts: ParaLies emits the shorter of two near-identical
    # copies, so which one it names is arbitrary.
    src_len = row.src_end - row.src_start
    dup_len = row.dup_end - row.dup_start
    src = covered(bed, row.src_chrom, row.src_start, row.src_end)
    dup = covered(bed, row.dup_chrom, row.dup_start, row.dup_end)
    return (src_len > 0 and src / src_len >= min_frac) or \
           (dup_len > 0 and dup / dup_len >= min_frac)


def match_synteny(truth, synteny_path):
    # A synteny row belongs to a block when one side is that block's uniquely
    # named contig and the other overlaps its source window.
    syn = pd.read_csv(synteny_path, sep="\t", dtype={"q_chr": str, "t_chr": str})
    by_contig = {}
    for r in syn.itertuples():
        for near, far, fs, fe in ((r.q_chr, r.t_chr, r.t_start, r.t_end),
                                  (r.t_chr, r.q_chr, r.q_start, r.q_end)):
            by_contig.setdefault(near, []).append((far, fs, fe, r.median_ks))

    out = []
    for b in truth.itertuples():
        best = None
        for far, fs, fe, ks in by_contig.get(b.dup_chrom, ()):
            if far != b.src_chrom or fe <= b.src_start or fs >= b.src_end:
                continue
            ov = min(fe, b.src_end) - max(fs, b.src_start)
            if best is None or ov > best[0]:
                best = (ov, ks)
        out.append(best[1] if best else float("nan"))
    return out
