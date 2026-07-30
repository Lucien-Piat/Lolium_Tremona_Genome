#!/bin/bash
# Usage : bash 01_remove_artefacts.sh

set -euo pipefail
SIF=images/sif/surgery.sif
REF=reference_data/lmultiflorum.tremona.fa
BED=results/dupclass/masked_intervals.bed
OUTDIR=clean_assembly
GAP=100
run() { singularity exec "${SIF}" "$@"; }
mkdir -p "${OUTDIR}"

PRIMARY="${OUTDIR}/lmultiflorum.tremona.primary.fa"
ALT="${OUTDIR}/lmultiflorum.tremona.alt.fa"
AGP="${OUTDIR}/lmultiflorum.tremona.primary.agp"
MAP="${OUTDIR}/coordinate_map.tsv"

[[ -s "${REF}.fai" ]] || run samtools faidx "${REF}"

# cut the broad intervals : kept sequence -> primary (gapped with 100 N), removed -> alt
run python3 scripts/04_surgery/01_remove_artefacts.py \
    --fasta "${REF}" \
    --bed "${BED}" \
    --gap "${GAP}" \
    --out-primary "${PRIMARY}" \
    --out-alt "${ALT}" \
    --out-agp "${AGP}" \
    --out-map "${MAP}"

run samtools faidx "${PRIMARY}"
run samtools faidx "${ALT}"
run samtools dict "${PRIMARY}" -o "${PRIMARY%.fa}.dict"

echo "primary -> ${PRIMARY}"
echo "alt     -> ${ALT}"
echo "agp     -> ${AGP}"
echo "map     -> ${MAP}"