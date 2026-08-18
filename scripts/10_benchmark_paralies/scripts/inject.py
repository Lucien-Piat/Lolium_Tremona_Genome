import math
import random
import re
from collections import defaultdict


_BASES = "ACGT"
_COMP  = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_ID_RE = re.compile(r"\b(ID|Parent)=([^;]+)")


from Bio.Data.CodonTable import standard_dna_table
_FWD   = standard_dna_table.forward_table
_STOPS = set(standard_dna_table.stop_codons)


def _aa(codon):
    if codon in _STOPS:
        return "*"
    return _FWD.get(codon)


# Observed divergence that a JC distance d produces: the inverse of the
# correction NG86 applies, so a sequence evolved at this rate is estimated
# back at exactly d.
def _jc_forward(d):
    return 0.75 * (1.0 - math.exp(-4.0 / 3.0 * d))


# Per codon position the three alternatives are split into synonymous and
# non-synonymous and applied at the rate NG86 inverts; a change that would
# create a premature stop is never offered.
def _evolve_cds(cds, ks, ka):
    p_s = _jc_forward(ks)
    p_n = _jc_forward(ka)
    out = list(cds)
    for c in range(len(cds) // 3):
        codon = cds[3 * c:3 * c + 3]
        ref_aa = _aa(codon)
        if ref_aa is None or ref_aa == "*":
            continue
        for i in range(3):
            syn, non = [], []
            for b in _BASES:
                if b == codon[i]:
                    continue
                alt_aa = _aa(codon[:i] + b + codon[i + 1:])
                if alt_aa is None or alt_aa == "*":
                    continue
                (syn if alt_aa == ref_aa else non).append(b)
            q_s = p_s * len(syn) / 3.0
            q_n = p_n * len(non) / 3.0
            r = random.random()
            if r < q_s:
                out[3 * c + i] = random.choice(syn)
            elif r < q_s + q_n:
                out[3 * c + i] = random.choice(non)
    return "".join(out)


def _indel_len():
    return max(1, int(random.expovariate(1.0 / INDEL_MEAN)) + 1)


P = snakemake.params

def _param(name, default):
    return type(default)(getattr(P, name, default))

SEED         = _param("seed", 42)
ARTEFACT_KS  = [float(x) for x in P.aks]
PARALOGUE_KS = [float(x) for x in P.pks]
GENE_COUNTS  = [int(x) for x in P.gene_counts]
N_REP        = _param("n", 10)
OMEGA        = _param("omega", 0.2)
INDEL_RATIO  = _param("indel_ratio", 0.12)
INDEL_MEAN   = _param("indel_mean", 3.0)
FLANK        = _param("flank", 2000)
BUFFER_GENES = _param("buffer_genes", 5)

random.seed(SEED)


sequences = {}
name = None
chunks = []
with open(snakemake.input.fa) as fh:
    for line in fh:
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
if name is not None:
    sequences[name] = "".join(chunks).upper()

total_bp = sum(len(s) for s in sequences.values())


gene_span = {}
mrna_gene = {}
cds_of    = defaultdict(list)
gff_rows  = defaultdict(list)

with open(snakemake.input.gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) != 9 or p[0] not in sequences:
            continue
        chrom = p[0]
        gff_rows[chrom].append(p)
        attrs = dict(_ID_RE.findall(p[8]))
        if p[2] == "gene":
            if "ID" in attrs:
                gene_span[(chrom, attrs["ID"])] = (int(p[3]), int(p[4]), p[6])
        elif p[2] == "mRNA":
            if "ID" in attrs and "Parent" in attrs:
                mrna_gene[(chrom, attrs["ID"])] = attrs["Parent"].split(",")[0]
        elif p[2] == "CDS":
            if "Parent" in attrs:
                cds_of[(chrom, attrs["Parent"].split(",")[0])].append(
                    (int(p[3]), int(p[4])))


best_mrna = {}
for (chrom, mrna), gene in mrna_gene.items():
    if (chrom, mrna) not in cds_of:
        continue
    span = sum(e - s + 1 for s, e in cds_of[(chrom, mrna)])
    key = (chrom, gene)
    if key not in best_mrna or span > best_mrna[key][1]:
        best_mrna[key] = (mrna, span)


usable = [
    (chrom, gid, gene_span[(chrom, gid)][0], gene_span[(chrom, gid)][1],
     gene_span[(chrom, gid)][2], best_mrna[(chrom, gid)][0])
    for (chrom, gid) in gene_span
    if (chrom, gid) in best_mrna and best_mrna[(chrom, gid)][1] >= 3 * 30
]


usable.sort(key=lambda g: (g[0], g[2]))
for chrom in sorted(sequences):
    n = sum(1 for g in usable if g[0] == chrom)


requests = [("artefact", ks, n) for ks in ARTEFACT_KS
            for n in GENE_COUNTS for _ in range(N_REP)]
requests += [("paralogue", ks, n) for ks in PARALOGUE_KS
             for n in GENE_COUNTS for _ in range(N_REP)]
random.shuffle(requests)

need = sum(n for _, _, n in requests) + BUFFER_GENES * len(requests)
if need > len(usable):
    raise RuntimeError(
        f"grid needs {need:,} anchor genes but only {len(usable):,} are usable. "
        f"Reduce n_blocks_per_cell, block_gene_counts, or the Ks grid, or add "
        f"more chromosomes."
    )


slack   = len(usable) - need
weights = [random.random() for _ in requests]
wsum    = sum(weights) or 1.0
extra   = [int(w / wsum * slack) for w in weights]

cursor, plan = 0, []
for (kind, ks, n), pad in zip(requests, extra):


    while cursor + n <= len(usable) and usable[cursor][0] != usable[cursor + n - 1][0]:
        cursor += 1
    if cursor + n > len(usable):
        raise RuntimeError(
            f"ran out of genes laying block {len(plan) + 1} of {len(requests)}; "
            f"chromosome boundaries consumed more slack than expected"
        )
    plan.append((kind, ks, n, cursor))
    cursor += n + BUFFER_GENES + pad


collapsed_contigs = []
hap_shared        = []
hap2_edits        = defaultdict(list)
injected_gff      = []
truth_rows        = []

for idx, (kind, ks, n_genes, gi) in enumerate(plan, start=1):
    block = usable[gi:gi + n_genes]
    chrom = block[0][0]
    seq = sequences[chrom]
    chrom_len = len(seq)

    g_start = min(g[2] for g in block)
    g_end   = max(g[3] for g in block)
    win_start = max(1, g_start - FLANK)
    win_end   = min(chrom_len, g_end + FLANK)
    s0, e0 = win_start - 1, win_end

    win = list(seq[s0:e0])
    in_cds = bytearray(len(win))


    for _c, gid, _gs, _ge, strand, mrna in block:
        segs = sorted(cds_of[(chrom, mrna)])
        pos = []
        if strand == "-":
            for s, e in sorted(segs, reverse=True):
                pos.extend(range(e, s - 1, -1))
        else:
            for s, e in segs:
                pos.extend(range(s, e + 1))
        if not pos or min(pos) < win_start or max(pos) > win_end:
            continue
        idxs = [p - win_start for p in pos]
        for k in idxs:
            in_cds[k] = 1
        cds = "".join(win[k] for k in idxs)
        if strand == "-":
            cds = cds.translate(_COMP)
        cds = cds[:len(cds) - len(cds) % 3]
        if len(cds) < 90:
            continue
        new = _evolve_cds(cds, ks, ks * OMEGA)
        for j, ch in enumerate(new):
            win[idxs[j]] = ch if strand == "+" else ch.translate(_COMP)


    p_sub   = _jc_forward(ks)
    p_indel = p_sub * INDEL_RATIO
    out, newpos, i = [], [0] * len(win), 0
    while i < len(win):
        newpos[i] = len(out)
        if in_cds[i]:
            out.append(win[i])
            i += 1
            continue
        r = random.random()
        if r < p_indel / 2:
            want, j = _indel_len(), i
            while j < len(win) and want > 0 and not in_cds[j]:
                newpos[j] = len(out)
                j += 1
                want -= 1
            i = j
            continue
        if r < p_indel:
            out.extend(random.choice(_BASES) for _ in range(_indel_len()))
        b = win[i]
        if b in _BASES and random.random() < p_sub:
            out.append(random.choice([x for x in _BASES if x != b]))
        else:
            out.append(b)
        i += 1

    dup_seq = "".join(out)
    name = f"{kind}_{idx:04d}_ks{ks}_g{n_genes}"


    for p in gff_rows[chrom]:
        fs, fe = int(p[3]), int(p[4])
        if fs < win_start or fe > win_end:
            continue
        q = list(p)
        q[0] = name
        q[3] = str(newpos[fs - win_start] + 1)
        q[4] = str(newpos[fe - win_start] + 1)
        if int(q[4]) < int(q[3]):
            continue
        q[8] = _ID_RE.sub(
            lambda m: f"{m.group(1)}=" + ",".join(f"{name}_{v}"
                                                  for v in m.group(2).split(",")),
            p[8])
        injected_gff.append("\t".join(q))


    collapsed_contigs.append((name, dup_seq))
    if kind == "artefact":
        hap2_edits[chrom].append((s0, e0, dup_seq))
    else:
        hap_shared.append((name, dup_seq))

    truth_rows.append(
        f"{name}\t{kind}\t{ks}\t{n_genes}\t{win_end - win_start + 1}\t"
        f"{chrom}\t{win_start}\t{win_end}\t{name}\t1\t{len(dup_seq)}"
    )


def _wrap(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


with open(snakemake.output.truth, "w") as fh:
    fh.write("region_id\ttype\tinjected_ks\tn_genes\tblock_bp\t"
             "src_chrom\tsrc_start\tsrc_end\tdup_chrom\tdup_start\tdup_end\n")
    fh.write("\n".join(truth_rows) + "\n")

with open(snakemake.output.collapsed, "w") as fh:
    for chrom in sorted(sequences):
        _wrap(fh, chrom, sequences[chrom])
    for name, seq in collapsed_contigs:
        _wrap(fh, name, seq)

with open(snakemake.output.collapsed_gff, "w") as fh:
    with open(snakemake.input.gff) as src:
        for line in src:
            fh.write(line)
    fh.write("\n".join(injected_gff) + "\n")
