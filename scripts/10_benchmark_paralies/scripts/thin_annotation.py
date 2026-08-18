import random
import re
from collections import defaultdict

drop = float(snakemake.wildcards.drop)
rng = random.Random(int(snakemake.params.seed) + int(drop * 1000))

ID_RE = re.compile(r"\bID=([^;]+)")
PARENT_RE = re.compile(r"\bParent=([^;]+)")

rows = []
gene_ids = []
with open(snakemake.input.gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) != 9:
            continue
        rows.append(p)
        if p[2] == "gene":
            m = ID_RE.search(p[8])
            if m:
                gene_ids.append(m.group(1))

n_drop = int(round(drop * len(gene_ids)))
doomed = set(rng.sample(gene_ids, n_drop)) if n_drop else set()


# Dropping a gene must take its mRNAs, exons and CDS with it.
children = defaultdict(list)
own_id = {}
for i, p in enumerate(rows):
    m = ID_RE.search(p[8])
    if m:
        own_id[i] = m.group(1)
    m = PARENT_RE.search(p[8])
    if m:
        for parent in m.group(1).split(","):
            children[parent].append(i)

removed_rows = set()
frontier = list(doomed)
while frontier:
    parent = frontier.pop()
    for i in children.get(parent, ()):
        if i in removed_rows:
            continue
        removed_rows.add(i)
        if i in own_id:
            frontier.append(own_id[i])
for i, p in enumerate(rows):
    if i in own_id and own_id[i] in doomed:
        removed_rows.add(i)

kept_genes = 0
with open(snakemake.output[0], "w") as out:
    for i, p in enumerate(rows):
        if i in removed_rows:
            continue
        out.write("\t".join(p) + "\n")
        if p[2] == "gene":
            kept_genes += 1

if kept_genes == 0:
    raise SystemExit("dropped every gene; ParaLies has nothing to run on")
