#!/bin/bash
# Usage : bash 02_metrics.sh
# Per-chromosome gene counts and basic metrics from the merged annotation.

set -euo pipefail
GFF=results/06_gene_annotation/tremona.merged.gff
FAI=reference_data/lmultiflorum.tremona.fa.fai
OUTDIR=results/06_gene_annotation
OUT="${OUTDIR}/annotation_metrics.tsv"

# per chromosome gene count
: > "${OUT}.body"
total=0
while read -r chr len _; do
    n=$(awk -F'\t' -v c="${chr}" '$1==c && $3=="gene"' "${GFF}" | wc -l)
    dens=$(awk -v n="${n}" -v l="${len}" 'BEGIN{printf "%.2f", (l>0)? n/(l/1e6):0}')
    printf "%s\t%d\t%d\t%s\n" "${chr}" "${len}" "${n}" "${dens}" >> "${OUT}.body"
    total=$(( total + n ))
done < "${FAI}"

# unplaced contigs
placed=$(awk -F'\t' '$2>=50000000{print $1}' "${FAI}" | sort -u)
total_all=$(awk -F'\t' '$3=="gene"' "${GFF}" | wc -l)
on_placed=$(awk -F'\t' '$3=="gene"{print $1}' "${GFF}" \
    | grep -Fxf <(echo "${placed}") | wc -l || true)
unplaced=$(( total_all - on_placed ))

{
    printf "chrom\tlength_bp\tgenes\tgenes_per_Mb\n"
    sort -k1,1 "${OUT}.body"
} > "${OUT}"
rm -f "${OUT}.body"
