#!/bin/bash
# Usage : bash 04_orientation_batch.sh [threads]

set -euo pipefail
SIF=images/sif/surgery.sif
ASM=clean_assembly/lmultiflorum.tremona.primary.fa
OUTDIR=results/orientation
THREADS=${1:-4}
SEG=500000          # large blocks only
MINID=85
MINBLOCK=300000
MINCHR=50000000
run() { singularity exec "${SIF}" "$@"; }
mkdir -p "${OUTDIR}"

REFS=(
    "Tremona_oriented    clean_assembly/lmultiflorum.tremona.primary.oriented.fa"
)

[[ -s "${ASM}.fai" ]] || run samtools faidx "${ASM}"

for entry in "${REFS[@]}"; do
    name=$(echo "${entry}" | awk '{print $1}')
    ref=$(echo "${entry}"  | awk '{print $2}')
    echo "=== ${name} : ${ref} ==="

    [[ -s "${ref}.fai" ]] || run samtools faidx "${ref}"

    run mashmap -r "${ref}" -q "${ASM}" \
        -s "${SEG}" --pi "${MINID}" -t "${THREADS}" \
        -o "${OUTDIR}/trem_vs_${name}.paf"

    run python3 scripts/04_surgery/02_orientation_dotplot.py \
        --paf "${OUTDIR}/trem_vs_${name}.paf" \
        --qfai "${ASM}.fai" \
        --rfai "${ref}.fai" \
        --min-chrom-len "${MINCHR}" \
        --min-block "${MINBLOCK}" \
        --title "Tremona vs ${name}" \
        --out-prefix "${OUTDIR}/orient_${name}"

done



