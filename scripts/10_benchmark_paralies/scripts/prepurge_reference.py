import bisect
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

fa_in = snakemake.input.fa
gff_in = snakemake.input.gff
bed_in = snakemake.input.bed
fa_out = snakemake.output.fa
gff_out = snakemake.output.gff

removed = common.read_bed(bed_in)


def wrap(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


sequences, order = {}, []
name, chunks = None, []
with open(fa_in) as fh:
    for line in fh:
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks)
            name = line[1:].split()[0]
            order.append(name)
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)


starts, cum_before, cut_bp = {}, {}, {}
for chrom, seq in sequences.items():
    ivs = removed.get(chrom, [])
    starts[chrom] = [s for s, _ in ivs]
    running, cum = 0, []
    for s, e in ivs:
        cum.append(running)
        running += e - s
    cum_before[chrom] = cum
    cut_bp[chrom] = running

with open(fa_out, "w") as out:
    for chrom in order:
        seq = sequences[chrom]
        ivs = removed.get(chrom, [])
        if not ivs:
            wrap(out, chrom, seq)
            continue
        kept, prev = [], 0
        for s, e in ivs:
            kept.append(seq[prev:s])
            prev = e
        kept.append(seq[prev:])
        wrap(out, chrom, "".join(kept))


# Excising interior sequence shifts every downstream coordinate, so the
# annotation is remapped, and a feature overlapping a cut is dropped whole
# rather than left frameshifted.
def inside_removal(chrom, lo, hi):
    ivs = removed.get(chrom, [])
    i = bisect.bisect_right(starts[chrom], lo) - 1
    for j in (i, i + 1):
        if 0 <= j < len(ivs):
            s, e = ivs[j]
            if lo < e and hi > s:
                return True
    return False


def remap(chrom, pos0):
    i = bisect.bisect_right(starts[chrom], pos0) - 1
    if i < 0:
        return pos0
    return pos0 - (cum_before[chrom][i] + (removed[chrom][i][1]
                                           - removed[chrom][i][0]))


ID_RE = re.compile(r"\bID=([^;]+)")
PARENT_RE = re.compile(r"\bParent=([^;]+)")

rows = []
with open(gff_in) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) == 9 and p[0] in sequences:
            rows.append(p)

doomed_ids = set()
for p in rows:
    lo, hi = int(p[3]) - 1, int(p[4])
    if inside_removal(p[0], lo, hi):
        m = ID_RE.search(p[8])
        if m:
            doomed_ids.add(m.group(1))


changed = True
while changed:
    changed = False
    for p in rows:
        m = ID_RE.search(p[8])
        mp = PARENT_RE.search(p[8])
        if not m or not mp:
            continue
        if m.group(1) in doomed_ids:
            continue
        if any(par in doomed_ids for par in mp.group(1).split(",")):
            doomed_ids.add(m.group(1))
            changed = True

kept_rows, n_genes, n_dropped = [], 0, 0
for p in rows:
    m = ID_RE.search(p[8])
    mp = PARENT_RE.search(p[8])
    gone = (m and m.group(1) in doomed_ids) or \
           (mp and any(par in doomed_ids for par in mp.group(1).split(",")))
    lo, hi = int(p[3]) - 1, int(p[4])
    if gone or inside_removal(p[0], lo, hi):
        if p[2] == "gene":
            n_dropped += 1
        continue
    p[3] = str(remap(p[0], lo) + 1)
    p[4] = str(remap(p[0], hi - 1) + 1)
    kept_rows.append(p)
    if p[2] == "gene":
        n_genes += 1

with open(gff_out, "w") as out:
    for p in kept_rows:
        out.write("\t".join(p) + "\n")
