import random
import re
from collections import defaultdict

ID_RE = re.compile(r"\bID=([^;]+)")
PARENT_RE = re.compile(r"\bParent=([^;]+)")

mode = snakemake.wildcards.mode
frac = float(snakemake.wildcards.frac)
random.seed(int(snakemake.params.seed))

rows, header = [], []
with open(snakemake.input.gff) as fh:
    for line in fh:
        if line.startswith("#"):
            header.append(line)
            continue
        if not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) == 9:
            rows.append(p)

gene_ids = [ID_RE.search(p[8]).group(1) for p in rows
            if p[2] == "gene" and ID_RE.search(p[8])]
n_target = int(round(len(gene_ids) * frac))
targets = set(random.sample(gene_ids, n_target)) if n_target else set()

if mode == "drop":


    children = defaultdict(list)
    id_of = {}
    for i, p in enumerate(rows):
        m = ID_RE.search(p[8])
        if m:
            id_of[m.group(1)] = i
        mp = PARENT_RE.search(p[8])
        if mp:
            for parent in mp.group(1).split(","):
                children[parent].append(i)

    doomed, stack = set(), list(targets)
    while stack:
        node = stack.pop()
        for i in children.get(node, ()):
            if i in doomed:
                continue
            doomed.add(i)
            m = ID_RE.search(rows[i][8])
            if m:
                stack.append(m.group(1))
    for gid in targets:
        if gid in id_of:
            doomed.add(id_of[gid])

    kept = [p for i, p in enumerate(rows) if i not in doomed]

# `frac` differs by mode: proportion of genes in drop, severity per gene in
# frag. Truncating a proportion of genes by a fixed amount costs no anchors
# at all -- a 1.2 kb gene stays far above the 90 bp floor.
elif mode == "frag":


    strand_of, cds_of_mrna = {}, defaultdict(list)
    for i, p in enumerate(rows):
        mp = PARENT_RE.search(p[8])
        mi = ID_RE.search(p[8])
        if p[2] == "mRNA" and mi:
            strand_of[mi.group(1)] = p[6]
        elif p[2] == "CDS" and mp:
            cds_of_mrna[mp.group(1).split(",")[0]].append(i)

    doomed = set()
    n_trimmed, n_emptied = 0, 0
    for mrna, idx in cds_of_mrna.items():
        minus = strand_of.get(mrna, "+") == "-"

        idx = sorted(idx, key=lambda i: int(rows[i][3]), reverse=minus)
        total = sum(int(rows[i][4]) - int(rows[i][3]) + 1 for i in idx)
        keep_bp = int(total * (1.0 - frac)) // 3 * 3
        if keep_bp >= total:
            continue
        if keep_bp < 3:
            doomed.update(idx)
            n_emptied += 1
            continue
        run = 0
        for i in idx:
            seg = int(rows[i][4]) - int(rows[i][3]) + 1
            if run >= keep_bp:
                doomed.add(i)
            elif run + seg > keep_bp:

                cut = keep_bp - run
                if minus:
                    rows[i][3] = str(int(rows[i][4]) - cut + 1)
                else:
                    rows[i][4] = str(int(rows[i][3]) + cut - 1)
            run += seg
        n_trimmed += 1

    kept = [p for i, p in enumerate(rows) if i not in doomed]

else:
    raise ValueError(f"unknown degradation mode {mode!r}")

with open(snakemake.output[0], "w") as fh:
    fh.writelines(header)
    for p in kept:
        fh.write("\t".join(p) + "\n")
