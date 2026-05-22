#!/usr/bin/env python3
"""Drop proteins shorter than 30 aa or with only X"""
import sys
from Bio import SeqIO

def main(in_faa, out_faa, min_len=30):
    kept = dropped = 0
    with open(out_faa, "w") as out:
        for rec in SeqIO.parse(in_faa, "fasta"):
            seq = str(rec.seq).replace(".", "").replace("*", "")
            if len(seq) < min_len or set(seq.upper()) <= {"X"}:
                dropped += 1
                continue
            rec.seq = type(rec.seq)(seq)
            SeqIO.write(rec, out, "fasta")
            kept += 1
    print(f"Kept {kept}, dropped {dropped}", file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])