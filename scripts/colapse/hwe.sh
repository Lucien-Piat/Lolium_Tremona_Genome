#!/bin/bash

# hwe + mosedepth

set -euo pipefail
SIF=images/sif/collapse_diag.sif
BAMDIR=mapping
OUTDIR=collapse_diag
SNPS="${OUTDIR}/cohort.snps.vcf.gz"
REGION=${1:-}
T=4
WIN=50000
run() { singularity exec "${SIF}" "$@"; }

run vcftools --gzvcf "${SNPS}" --hardy --out "${OUTDIR}/trem"

CHROM_ARG=""; [[ -n "${REGION}" ]] && CHROM_ARG="--chrom ${REGION%%:*}"
DEPTHLIST="${OUTDIR}/depth_files.txt"; : > "${DEPTHLIST}"
shopt -s nullglob
bams=("${BAMDIR}"/TREM*.dedup.bam)
if (( ${#bams[@]} == 0 )); then
    echo "ATTENTION : aucun BAM dans ${BAMDIR}/, figures en FIS seul."
else
    for b in "${bams[@]}"; do
        s=$(basename "${b}" .dedup.bam)
        [[ -s "${b}.bai" ]] || run samtools index -@ "${T}" "${b}"
        run mosdepth -t "${T}" -n --fast-mode ${CHROM_ARG} -b "${WIN}" "${OUTDIR}/${s}" "${b}"
        echo "${s} ${OUTDIR}/${s}.regions.bed.gz" >> "${DEPTHLIST}"
    done
fi

echo "Termine -> ${OUTDIR}/"