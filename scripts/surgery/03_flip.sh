#!/bin/bash
# Usage : bash 05_flip_chromosomes.sh chr1,chr4,chr5
# Reverse-complement the listed chromosomes in the primary, keep order and names.

set -euo pipefail
SIF=images/sif/surgery.sif
ASM=clean_assembly/lmultiflorum.tremona.primary.fa
OUT=clean_assembly/lmultiflorum.tremona.primary.oriented.fa
FLIP=${1:?give a comma-separated flip list, e.g. chr1,chr4,chr5}
run() { singularity exec "${SIF}" "$@"; }

FLIP_SET=" ${FLIP//,/ } "

[[ -s "${ASM}.fai" ]] || run samtools faidx "${ASM}"

: > "${OUT}"
while read -r chr _; do
    if [[ "${FLIP_SET}" == *" ${chr} "* ]]; then
        echo "flipping ${chr}"
        run samtools faidx "${ASM}" "${chr}" | run seqkit seq -r -p >> "${OUT}"
    else
        run samtools faidx "${ASM}" "${chr}" >> "${OUT}"
    fi
done < "${ASM}.fai"

run samtools faidx "${OUT}"
echo "oriented fasta -> ${OUT}"