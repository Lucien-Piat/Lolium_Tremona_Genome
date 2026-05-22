
#!/usr/bin/env python3
"""Convert Oatk FASTA + BED to GenBank format."""
import sys
from Bio import SeqIO # type: ignore
from Bio.SeqFeature import SeqFeature, FeatureLocation # type: ignore

def feature_type(gene_name):
    if gene_name.startswith('trn'):
        return 'tRNA'
    if gene_name.startswith('rrn'):
        return 'rRNA'
    return 'CDS'

def bed_to_features(bed_path):
    features = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip().split('\t')
            seq, start, end, gene, _, strand = parts[:6]
            feat = SeqFeature(
                FeatureLocation(int(start), int(end),
                                strand=1 if strand == '+' else -1),
                type=feature_type(gene),
                qualifiers={'gene': [gene], 'product': [gene]},
            )
            features.setdefault(seq, []).append(feat)
    return features

def convert(fasta_path, bed_path, out_path, organism, organelle):
    features = bed_to_features(bed_path)
    records = []
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        rec.annotations['molecule_type'] = 'DNA'
        rec.annotations['organism'] = organism
        rec.annotations['topology'] = 'linear'
        rec.description = f'{organism} {organelle}, contig {rec.id}'
        rec.features.extend(features.get(rec.id, []))
        records.append(rec)
    SeqIO.write(records, out_path, 'genbank')
    print(f'Wrote {len(records)} record(s) to {out_path}')

if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit('usage: bed_to_gb.py FASTA BED OUT.gb ORGANISM ORGANELLE')
    convert(*sys.argv[1:6])