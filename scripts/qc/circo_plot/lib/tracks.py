"""Data preparation for the Circos plot tracks."""
from bisect import bisect_left
from collections import OrderedDict, defaultdict

from Bio import SeqIO # type: ignore


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
    print("Finding gaps")
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
    print("Computing G/C")
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
    print("Computing gene density")
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
    print("Computing busco list")
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
    """Lenient GenBank feature parser.

    Skips individual malformed features (not whole records). For locations
    like 'join(134971..135351,1..500)' or even truncated 'join(134971..135351'
    we just take the min/max of any numeric ranges found.

    Returns [(chrom_id, name, start, end, strand, ftype), ...].
    """
    print("Computing organellear features")
    import re

    out = []
    chrom_id = None
    in_features = False
    ftype, floc, fquals = None, None, {}

    header_re = re.compile(r"^ {5}(\S+)\s+(.+?)\s*$")

    def flush():
        nonlocal ftype, floc, fquals
        if ftype in types and floc:
            parsed = _parse_loc(floc)
            if parsed is not None:
                start, end, strand = parsed
                name = fquals.get("gene") or fquals.get("product") or "?"
                out.append((chrom_id, name, start, end, strand, ftype))
        ftype, floc, fquals = None, None, {}

    with open(gb_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("LOCUS"):
                flush()
                parts = line.split()
                chrom_id = parts[1] if len(parts) > 1 else None
                in_features = False
            elif line.startswith("FEATURES"):
                in_features = True
            elif line.startswith("ORIGIN") or line.startswith("//"):
                flush()
                in_features = False
            elif in_features:
                # New feature: exactly 5 leading spaces
                if len(line) - len(line.lstrip(" ")) == 5:
                    m = header_re.match(line)
                    if m:
                        flush()
                        ftype = m.group(1)
                        floc = m.group(2).strip()
                elif line.startswith(" " * 21):
                    content = line[21:].strip()
                    if content.startswith("/"):
                        qm = re.match(r"/(\w+)=?\"?([^\"]*)\"?", content)
                        if qm:
                            fquals[qm.group(1)] = qm.group(2)
                    elif floc is not None:
                        floc += content
    flush()
    return out


def _parse_loc(loc_str):
    """Parse a GenBank location string. Returns (start_0based, end, strand) or None."""
    import re
    strand = 1
    if loc_str.startswith("complement("):
        strand = -1
        loc_str = loc_str[11:]
    # Strip wrappers
    for prefix in ("join(", "order("):
        if loc_str.startswith(prefix):
            loc_str = loc_str[len(prefix):]
    matches = re.findall(r"(\d+)\.\.(\d+)", loc_str)
    if not matches:
        return None
    starts = [int(m[0]) for m in matches]
    ends   = [int(m[1]) for m in matches]
    s, e = min(starts) - 1, max(ends)
    return (s, e, strand) if e > s else None


# Links

def load_self_synteny(path):
    print("Loading chr syntheny")
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
    print("Loading organellar syntheny")
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

def organelle_gc(fasta_path, offsets, scale, window=2000):
    """GC content for organelle contigs, in virtual sector coordinates.
    Returns [(start_virt, end_virt, gc_pct), ...]."""
    print("Computing organellar G/C")
    out = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if rec.id not in offsets:
            continue
        off = offsets[rec.id]
        seq = str(rec.seq).upper()
        for start in range(0, len(seq), window):
            end = min(start + window, len(seq))
            chunk = seq[start:end]
            valid = sum(1 for b in chunk if b in "ACGT")
            if valid == 0:
                continue
            gc = sum(1 for b in chunk if b in "GC") / valid * 100
            out.append(((off + start) * scale, (off + end) * scale, gc))
    return out


def load_and_average_coverage(file_paths):
    """Reads multiple coverage files and averages the values per bin.
    Returns [(chrom, start, end, avg_coverage), ...]."""
    print("Computing coverage bins")
    if not file_paths:
        return []

    sums = defaultdict(float)
    counts = defaultdict(int)

    for path in file_paths:
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    p = line.rstrip().split("\t")
                    if len(p) >= 4:
                        chrom = p[0]
                        start = int(p[1])
                        end = int(p[2])
                        val = float(p[3])
                        sums[(chrom, start, end)] += val
                        counts[(chrom, start, end)] += 1

    out = []
    for (chrom, start, end), total in sums.items():
        out.append((chrom, start, end, total / counts[(chrom, start, end)]))
    
    return out