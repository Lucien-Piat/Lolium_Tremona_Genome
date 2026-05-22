"""Data preparation for the Circos plot tracks."""
from bisect import bisect_left
from collections import OrderedDict, defaultdict

from Bio import SeqIO


# Karyotype

def read_fai(fai_path):
    """Return OrderedDict of {chrom: length} from a .fai file."""
    chroms = OrderedDict()
    with open(fai_path) as fh:
        for line in fh:
            name, length, *_ = line.rstrip().split("\t")
            chroms[name] = int(length)
    return chroms


def organelle_lengths(*fasta_paths):
    out = OrderedDict()
    for path in fasta_paths:
        for rec in SeqIO.parse(path, "fasta"):
            out[rec.id] = len(rec.seq)
    return out


def find_gaps(fasta_path, min_gap=1000):
    """Find runs of N of at least min_gap bp. Returns [(chrom, start, end), ...]."""
    gaps = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper()
        in_gap, gstart = False, 0
        for i, base in enumerate(seq):
            if base == "N":
                if not in_gap:
                    in_gap, gstart = True, i
            elif in_gap:
                if i - gstart >= min_gap:
                    gaps.append((rec.id, gstart, i))
                in_gap = False
        if in_gap and len(seq) - gstart >= min_gap:
            gaps.append((rec.id, gstart, len(seq)))
    return gaps


# GC content

def gc_windows(fasta_path, window=1_000_000):
    """Return [(chrom, start, end, gc_pct), ...]. Ns excluded from denominator."""
    out = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper()
        for start in range(0, len(seq), window):
            end = min(start + window, len(seq))
            chunk = seq[start:end]
            valid = sum(1 for b in chunk if b in "ACGT")
            if valid == 0:
                continue
            gc = sum(1 for b in chunk if b in "GC") / valid * 100
            out.append((rec.id, start, end, gc))
    return out


# Gene density

def gene_density(gff_path, chrom_lengths, window=1_000_000):
    """Return [(chrom, start, end, n_genes), ...]."""
    starts = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) >= 9 and parts[2] == "gene":
                starts[parts[0]].append(int(parts[3]))
    out = []
    for chrom, length in chrom_lengths.items():
        ss = sorted(starts.get(chrom, []))
        for ws in range(0, length, window):
            we = min(ws + window, length)
            lo = bisect_left(ss, ws)
            hi = bisect_left(ss, we)
            out.append((chrom, ws, we, hi - lo))
    return out


# BUSCO orthologs

def busco_orthologs(table_path):
    """Parse BUSCO full table. Returns [(chrom, midpoint, status), ...] excluding Missing."""
    out = []
    with open(table_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip().split("\t")
            if len(p) < 5 or p[1] == "Missing":
                continue
            try:
                pos = (int(p[3]) + int(p[4])) // 2
            except ValueError:
                continue
            out.append((p[2], pos, p[1]))
    return out


# Organelle annotations

def organelle_features(gb_path, types=("CDS", "tRNA", "rRNA")):
    """Parse GenBank, return [(chrom, name, start, end, strand, type), ...]."""
    out = []
    for rec in SeqIO.parse(gb_path, "genbank"):
        for feat in rec.features:
            if feat.type not in types:
                continue
            name = feat.qualifiers.get("gene",
                   feat.qualifiers.get("product", ["?"]))[0]
            out.append((
                rec.id, name,
                int(feat.location.start), int(feat.location.end),
                feat.location.strand or 1, feat.type,
            ))
    return out


# Links

def load_self_synteny(path):
    out = []
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) < 7:
                continue
            out.append({
                "q_chr": p[0], "q_start": int(p[1]), "q_end": int(p[2]),
                "t_chr": p[3], "t_start": int(p[4]), "t_end": int(p[5]),
                "strand": p[6],
                "n_genes": int(p[7]) if len(p) > 7 else 1,
            })
    return out


def load_organelle_links(path):
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) < 9:
                continue
            out.append({
                "org_chr": p[0], "org_start": int(p[1]), "org_end": int(p[2]),
                "nuc_chr": p[3], "nuc_start": int(p[4]), "nuc_end": int(p[5]),
                "strand": p[6], "tag": p[7], "identity": float(p[8]),
            })
    return out