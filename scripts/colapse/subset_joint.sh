#!/bin/bash
# Usage : bash 01_joint_call.sh

set -euo pipefail
SIF=images/sif/collapse_diag.sif
VCFDIR=vcf
REF=reference_data/lmultiflorum.tremona.fa
OUTDIR=collapse_diag
REGION=${1:-}
run() { singularity exec "${SIF}" "$@"; }
mkdir -p "${OUTDIR}"

[[ -s "${REF}.fai" ]] || run samtools faidx "${REF}"
DICT="${REF%.fa}.dict"
[[ -s "${DICT}" ]] || run gatk CreateSequenceDictionary -R "${REF}" -O "${DICT}"

L_ARG=""; [[ -n "${REGION}" ]] && L_ARG="-L ${REGION}"

V_ARGS=""; for g in "${VCFDIR}"/TREM*.g.vcf.gz; do V_ARGS+=" -V ${g}"; done
run gatk --java-options "-Xmx8G" CombineGVCFs -R "${REF}" ${V_ARGS} ${L_ARG} \
    -O "${OUTDIR}/cohort.g.vcf.gz"

run gatk --java-options "-Xmx8G" GenotypeGVCFs -R "${REF}" \
    -V "${OUTDIR}/cohort.g.vcf.gz" ${L_ARG} -O "${OUTDIR}/cohort.vcf.gz"

run bash -c "bcftools view -v snps -m2 -M2 '${OUTDIR}/cohort.vcf.gz' \
    -Oz -o '${OUTDIR}/cohort.snps.vcf.gz'"
run bcftools index -t "${OUTDIR}/cohort.snps.vcf.gz"

rm -f "${OUTDIR}/cohort.g.vcf.gz" "${OUTDIR}/cohort.g.vcf.gz.tbi"

echo "VCF cohorte pret -> ${OUTDIR}/cohort.snps.vcf.gz"