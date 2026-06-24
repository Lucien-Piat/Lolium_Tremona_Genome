"""Data preparation for the Circos plot tracks."""
import gzip
import re
from bisect import bisect_left
from collections import OrderedDict, defaultdict

import numpy as np  # type: ignore
from Bio import SeqIO # type: ignore

# Karyotype
def read_fai(fai_path):
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
    print("Computing organellear features")
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
    strand = 1
    if loc_str.startswith("complement("):
        strand = -1
        loc_str = loc_str[11:]
    for prefix in ("join(", "order("):
        if loc_str.startswith(prefix):
            loc_str = loc_str[len(prefix):]
    matches = re.findall(r"(\d+)\.\.(\d+)", loc_str)
    if not matches: return None
    starts = [int(m[0]) for m in matches]
    ends   = [int(m[1]) for m in matches]
    s, e = min(starts) - 1, max(ends)
    return (s, e, strand) if e > s else None

# Links
def load_self_synteny(path):
    print("Loading chr syntheny")
    out = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) < 7: continue
            out.append({
                "q_chr": p[0], "q_start": int(p[1]), "q_end": int(p[2]),
                "t_chr": p[3], "t_start": int(p[4]), "t_end": int(p[5]),
                "strand": p[6], "n_genes": int(p[7]) if len(p) > 7 else 1,
            })
    return out

def load_organelle_links(path):
    print("Loading organellar syntheny")
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip().split("\t")
            if len(p) < 9: continue
            out.append({
                "org_chr": p[0], "org_start": int(p[1]), "org_end": int(p[2]),
                "nuc_chr": p[3], "nuc_start": int(p[4]), "nuc_end": int(p[5]),
                "strand": p[6], "tag": p[7], "identity": float(p[8]),
            })
    return out

def organelle_gc(fasta_path, offsets, scale, window=2000):
    out = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if rec.id not in offsets: continue
        off = offsets[rec.id]
        seq = str(rec.seq).upper()
        for start in range(0, len(seq), window):
            end = min(start + window, len(seq))
            chunk = seq[start:end]
            valid = sum(1 for b in chunk if b in "ACGT")
            if valid == 0: continue
            gc = sum(1 for b in chunk if b in "GC") / valid * 100
            out.append(((off + start) * scale, (off + end) * scale, gc))
    return out

# TE Mapping and Categorization
def load_te_mapping(path):
    print("Loading TE mapping table")
    mapping = {}
    with open(path) as fh:
        next(fh) # Skip header
        for line in fh:
            p = line.strip().split('\t')
            if len(p) >= 3:
                family_id = p[1]
                te_class = p[2]
                mapping[family_id] = te_class
    return mapping

def categorize_te(attr_str, mapping_dict):
    """Map Motif to precise color categories."""
    m = re.search(r'Motif:([^"\s]+)', attr_str)
    te_class = None
    if m:
        clean_motif = m.group(1).replace("-int", "").replace("_I", "")
        if clean_motif in mapping_dict:
            te_class = mapping_dict[clean_motif]
            
    s = (te_class if te_class else attr_str).upper()
    
    if "GYPSY" in s: return "Gypsy"
    if "COPIA" in s: return "Copia"
    if "PAO" in s or "ERV" in s or "LTR" in s: return "other LTR"
    if "L1" in s or "RTE" in s or "LINE" in s or "NON_LTR" in s or "NON-LTR" in s: return "LINE"
    if "SINE" in s or "TRNA" in s: return "SINE"
    if "HELITRON" in s or "TCMAR" in s or "MULE" in s or "HAT" in s or "HARBINGER" in s or "CMC" in s or "TIR" in s or "DNA" in s: return "DNA transposon"
    if "UNKNOWN" in s: return "Unknown"
    
    return "Other"

# TE Density
def te_density(gff_path, chrom_lengths, mapping_dict, window=1_000_000):
    print("Computing TE density bins")
    cov = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.rstrip('\n').split("\t")
            if len(p) < 9: continue
            chrom, start, end = p[0], int(p[3]), int(p[4])
            if chrom not in chrom_lengths: continue

            fam = categorize_te(p[8], mapping_dict)
            win_idx = start // window
            cov[chrom][win_idx][fam] += (end - start)

    out = {}
    for chrom, length in chrom_lengths.items():
        n_bins = (length + window - 1) // window 
        x = []
        widths = []
        for i in range(n_bins):
            ws = i * window
            we = min(ws + window, length)
            x.append((ws + we) / 2)           
            widths.append((we - ws) * 0.85)   
            
        x = np.array(x)
        widths = np.array(widths)
        
        fams_in_chrom = set()
        for w in range(n_bins):
            fams_in_chrom.update(cov[chrom][w].keys())
            
        y_stack = {f: np.zeros(n_bins) for f in fams_in_chrom}
        for w in range(n_bins):
            ws = w * window
            we = min(ws + window, length)
            actual_width = we - ws
            if actual_width <= 0: continue
            
            for f, val in cov[chrom][w].items():
                y_stack[f][w] = (val / actual_width) * 100
                
        out[chrom] = (x, widths, y_stack)
        
    return out

# SNP Density
def snp_density(vcf_path, chrom_lengths, window=1_000_000):
    print("Computing SNP density")
    starts = defaultdict(list)
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.split("\t", 2)
            starts[p[0]].append(int(p[1]))
            
    out = []
    for chrom, length in chrom_lengths.items():
        ss = sorted(starts.get(chrom, []))
        for ws in range(0, length, window):
            we = min(ws + window, length)
            lo = bisect_left(ss, ws)
            hi = bisect_left(ss, we)
            out.append((chrom, ws, we, hi - lo))
    return out