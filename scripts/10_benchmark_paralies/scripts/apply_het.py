import math
import os
import random
import re
from collections import defaultdict

from Bio.Data.CodonTable import standard_dna_table

het = float(snakemake.wildcards.het)
p_site = het / 2.0
gff_path = getattr(snakemake.input, "gff", None)
seed_salt = int(getattr(snakemake.params, "salt", 0))

BASES = "ACGT"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = set(standard_dna_table.stop_codons)
ID_RE = re.compile(r"\b(ID|Parent)=([^;]+)")


rng = random.Random(int(snakemake.params.seed) * 1000 + seed_salt
                    + int(het * 1e6))


def wrap(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def load_cds(path, wanted):
    mrna_gene, segs = {}, defaultdict(list)
    strand = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[0] not in wanted:
                continue
            attrs = dict(ID_RE.findall(p[8]))
            if p[2] == "mRNA" and "ID" in attrs and "Parent" in attrs:
                mrna_gene[(p[0], attrs["ID"])] = attrs["Parent"].split(",")[0]
                strand[(p[0], attrs["ID"])] = p[6]
            elif p[2] == "CDS" and "Parent" in attrs:
                segs[(p[0], attrs["Parent"].split(",")[0])].append(
                    (int(p[3]), int(p[4])))
    best = {}
    for key, gene in mrna_gene.items():
        if key not in segs:
            continue
        span = sum(e - s + 1 for s, e in segs[key])
        gkey = (key[0], gene)
        if gkey not in best or span > best[gkey][1]:
            best[gkey] = (key, span)
    out = defaultdict(list)
    for (chrom, _gene), (key, _span) in best.items():
        out[chrom].append((key[1], strand[key], sorted(segs[key])))
    return out


def ordered_positions(strand, segs):
    pos = []
    if strand == "-":
        for s, e in sorted(segs, reverse=True):
            pos.extend(range(e, s - 1, -1))
    else:
        for s, e in segs:
            pos.extend(range(s, e + 1))
    return pos


# Substitutions only, and never one creating a stop: coordinates survive so
# a single annotation serves every het level.
def mutate_cds(cds, p):
    if p <= 0:
        return cds, 0
    out = list(cds)
    n_codons = len(cds) // 3
    changed = 0
    log1mp = math.log1p(-p)
    i, n = 0, n_codons * 3
    while i < n:
        u = rng.random()
        if u <= 0.0:
            break
        i += int(math.log(u) / log1mp)
        if i >= n:
            break
        c, o = divmod(i, 3)
        codon = "".join(out[3 * c:3 * c + 3])
        if codon in STOPS:
            i += 1
            continue
        alts = []
        for b in BASES:
            if b == codon[o]:
                continue
            if codon[:o] + b + codon[o + 1:] in STOPS:
                continue
            alts.append(b)
        if alts:
            out[i] = rng.choice(alts)
            changed += 1
        i += 1
    return "".join(out), changed


def mutate_free(seq_list, mask, p):
    if p <= 0:
        return 0
    changed = 0
    log1mp = math.log1p(-p)
    i, n = 0, len(seq_list)
    while i < n:
        u = rng.random()
        if u <= 0.0:
            break
        i += int(math.log(u) / log1mp)
        if i >= n:
            break
        if not mask[i]:
            b = seq_list[i]
            if b in BASES:
                seq_list[i] = rng.choice([x for x in BASES if x != b])
                changed += 1
        i += 1
    return changed


sequences, order = {}, []
name, chunks = None, []
with open(snakemake.input.fa if hasattr(snakemake.input, "fa")
          else snakemake.input[0]) as fh:
    for line in fh:
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            order.append(name)
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks).upper()

cds_by_chrom = load_cds(gff_path, set(sequences)) if gff_path else {}

n_cds, n_free, n_total = 0, 0, 0
with open(snakemake.output[0], "w") as out:
    for chrom in order:
        seq = list(sequences[chrom])
        n_total += len(seq)
        mask = bytearray(len(seq))

        for mrna, strand, segs in cds_by_chrom.get(chrom, ()):
            pos = ordered_positions(strand, segs)
            idxs = [p - 1 for p in pos]
            if not idxs or idxs[0] >= len(seq) or max(idxs) >= len(seq):
                continue
            for k in idxs:
                mask[k] = 1
            cds = "".join(seq[k] for k in idxs)
            if strand == "-":
                cds = cds.translate(COMP)
            usable = len(cds) - len(cds) % 3
            if usable < 3:
                continue
            new, changed = mutate_cds(cds[:usable], p_site)
            n_cds += changed
            for j, ch in enumerate(new):
                seq[idxs[j]] = ch if strand == "+" else ch.translate(COMP)

        n_free += mutate_free(seq, mask, p_site)
        wrap(out, chrom, "".join(seq))


# gffread trusts any .fai it finds beside the FASTA, so a stale one from an
# earlier run must not survive.
stale = snakemake.output[0] + ".fai"
if os.path.exists(stale):
    os.remove(stale)

realised = (n_cds + n_free) / n_total if n_total else 0.0
if het > 0 and n_total and abs(realised - p_site) > 0.25 * p_site:
    raise SystemExit(f"realised divergence {realised:.5f} != requested {p_site:.5f}")
