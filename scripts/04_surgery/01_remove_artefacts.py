#!/usr/bin/env python3

import argparse 
import pysam # type: ignore


def load_bed(path):
    regions = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if end > start:
                regions.setdefault(chrom, []).append((start, end))
    for chrom, ivs in regions.items():
        ivs.sort()
        merged = [ivs[0]]
        for s, e in ivs[1:]:
            ls, le = merged[-1]
            if s <= le:
                merged[-1] = (ls, max(le, e))
            else:
                merged.append((s, e))
        regions[chrom] = merged
    return regions


def complement(removals, length):
    keeps, pos = [], 0
    for s, e in removals:
        if s > pos:
            keeps.append((pos, s))
        pos = max(pos, e)
    if pos < length:
        keeps.append((pos, length))
    return keeps


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--gap", type=int, default=100)
    ap.add_argument("--out-primary", required=True)
    ap.add_argument("--out-alt", required=True)
    ap.add_argument("--out-agp", required=True)
    ap.add_argument("--out-map", required=True)
    args = ap.parse_args()

    fa = pysam.FastaFile(args.fasta)
    removals = load_bed(args.bed)
    gapseq = "N" * args.gap

    fp = open(args.out_primary, "w")
    fa_alt = open(args.out_alt, "w")
    agp = open(args.out_agp, "w")
    cmap = open(args.out_map, "w")
    cmap.write("chrom\told_start\told_end\tnew_start\tnew_end\n")

    for chrom in fa.references:
        length = fa.get_reference_length(chrom)
        rem = removals.get(chrom, [])

        for s, e in rem:
            fa_alt.write(">%s__%d_%d\n%s\n" % (chrom, s, e, wrap(fa.fetch(chrom, s, e))))

        keeps = complement(rem, length)
        if not keeps:
            continue 

        parts, new_pos, part_no = [], 0, 0
        for i, (s, e) in enumerate(keeps):
            if i > 0:
                parts.append(gapseq)
                part_no += 1
                agp.write("%s\t%d\t%d\t%d\tU\t%d\tscaffold\tno\tna\n" % (
                    chrom, new_pos + 1, new_pos + args.gap, part_no, args.gap))
                new_pos += args.gap
            seg_len = e - s
            parts.append(fa.fetch(chrom, s, e))
            part_no += 1
            agp.write("%s\t%d\t%d\t%d\tW\t%s\t%d\t%d\t+\n" % (
                chrom, new_pos + 1, new_pos + seg_len, part_no, chrom, s + 1, e))
            cmap.write("%s\t%d\t%d\t%d\t%d\n" % (chrom, s, e, new_pos, new_pos + seg_len))
            new_pos += seg_len

        fp.write(">%s\n%s\n" % (chrom, wrap("".join(parts))))

    fp.close()
    fa_alt.close()
    agp.close()
    cmap.close()


if __name__ == "__main__":
    main()