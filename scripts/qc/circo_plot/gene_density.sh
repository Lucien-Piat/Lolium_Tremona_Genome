#!/bin/bash

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/annotation.sif)
TARGET=$(readlink -f reference_data/lmultiflorum.tremona.full.fa)
PLACED=$(readlink -f reference_data/lmultiflorum.tremona.placed.fa)
BRUN_FA=$(readlink -f reference_data/LOLMU.fa)
BRUN_GFF=$(readlink -f reference_data/LOLMU.genes.matched.gff)
KYUSS_FA=$(readlink -f reference_data/kyuss_v2.fasta)
KYUSS_GFF=$(readlink -f reference_data/kyuss_v2.gff)
OUTDIR="results/07_annotation"
T=4
ROOT=$(pwd)
run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Brunharo
#echo "[$(date)] Liftoff pass 1: Brunharo"
#run liftoff \
#    -g "${BRUN_GFF}" \
#    -o tremona.brunharo.gff \
#    -u tremona.brunharo.unmapped.txt \
#    -dir liftoff_brunharo_tmp \
#    -p "${T}" \
#    -copies \
#    -polish \
#    -mm2_options="-a --end-bonus 5 --eqx -N 50 -p 0.5 -I 4G" \
#    "${TARGET}" "${BRUN_FA}"
#mv tremona.brunharo.gff_polished tremona.brunharo.gff
#rm -rf liftoff_brunharo_tmp

# Kyuss v2
echo "[$(date)] Liftoff pass 2: Kyuss v2"
run liftoff \
    -g "${KYUSS_GFF}" \
    -o tremona.kyuss.gff \
    -u tremona.kyuss.unmapped.txt \
    -dir liftoff_kyuss_tmp \
    -p "${T}" \
    -copies \
    -polish \
    -mm2_options="-a --end-bonus 5 --eqx -N 50 -p 0.5 -I 4G" \
    "${TARGET}" "${KYUSS_FA}"
mv tremona.kyuss.gff_polished tremona.kyuss.gff
rm -rf liftoff_kyuss_tmp

# Extract as BED 
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{print $1, $4-1, $5, $9}' \
    tremona.brunharo.gff > brunharo.genes.bed
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{print $1, $4-1, $5, $9}' \
    tremona.kyuss.gff > kyuss.genes.bed

# Find Kyuss genes with NO overlap to Brunharo 
run bedtools intersect -v -a kyuss.genes.bed -b brunharo.genes.bed \
    | awk -F'\t' '{split($4, a, ";"); for (i in a) if (a[i] ~ /^ID=/) {sub(/^ID=/, "", a[i]); print a[i]}}' \
    > kyuss.unique_gene_ids.txt

# Pull those genes (and their children) from the Kyuss GFF
run python3 -c "
import sys
keep = set(open('kyuss.unique_gene_ids.txt').read().split())
keep_lines = []
current_gene = None
for line in open('tremona.kyuss.gff'):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9:
        continue
    attrs = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
    if f[2] == 'gene':
        current_gene = attrs.get('ID')
        if current_gene in keep:
            keep_lines.append(line)
    else:
        parent = attrs.get('Parent', '').split(',')[0]
        if current_gene in keep:
            keep_lines.append(line)
open('kyuss.unique.gff', 'w').writelines(keep_lines)
"

# Concatenate
{
    echo "##gff-version 3"
    awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {
        if ($9 !~ /source_annot=/) $9 = $9";source_annot=brunharo"
        print
    }' tremona.brunharo.gff
    awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {
        if ($9 !~ /source_annot=/) $9 = $9";source_annot=kyuss"
        print
    }' kyuss.unique.gff
} > tremona.merged.gff

# placed-only version 
run seqkit seq -n -i "${PLACED}" > placed.chroms.txt
awk -F'\t' 'NR==FNR{a[$1]=1; next} /^#/ || a[$1]' \
    placed.chroms.txt tremona.merged.gff > tremona.merged.placed.gff

# Summary 
echo "[$(date)] Summary:" | tee mapping_stats.txt
{
    N_BRUN_REF=$(awk -F'\t' '$3=="gene"' "${BRUN_GFF}" | wc -l)
    N_BRUN_MAP=$(awk -F'\t' '$3=="gene"' tremona.brunharo.gff | wc -l)
    N_KYUSS_REF=$(awk -F'\t' '$3=="gene"' "${KYUSS_GFF}" | wc -l)
    N_KYUSS_MAP=$(awk -F'\t' '$3=="gene"' tremona.kyuss.gff | wc -l)
    N_KYUSS_KEPT=$(awk -F'\t' '$3=="gene"' kyuss.unique.gff | wc -l)
    N_MERGED=$(awk -F'\t' '$3=="gene"' tremona.merged.gff | wc -l)
    N_PLACED=$(awk -F'\t' '$3=="gene"' tremona.merged.placed.gff | wc -l)
    echo "Brunharo: ${N_BRUN_MAP}/${N_BRUN_REF} lifted ($(awk "BEGIN{printf \"%.1f\", ${N_BRUN_MAP}/${N_BRUN_REF}*100}")%)"
    echo "Kyuss:    ${N_KYUSS_MAP}/${N_KYUSS_REF} lifted ($(awk "BEGIN{printf \"%.1f\", ${N_KYUSS_MAP}/${N_KYUSS_REF}*100}")%)"
    echo "Kyuss rescuing in regions Brunharo missed: ${N_KYUSS_KEPT}"
    echo "Merged total: ${N_MERGED} genes"
    echo "On 7 placed chromosomes: ${N_PLACED} genes ($(awk "BEGIN{printf \"%.1f\", ${N_PLACED}/${N_MERGED}*100}")%)"
} | tee -a mapping_stats.txt


run pigz -p "${T}" tremona.brunharo.gff tremona.kyuss.gff kyuss.unique.gff
rm -f brunharo.genes.bed kyuss.genes.bed kyuss.unique_gene_ids.txt placed.chroms.txt
