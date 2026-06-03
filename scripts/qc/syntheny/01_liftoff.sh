#!/bin/bash

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/annotation.sif)
REF_FA=$(readlink -f reference_data/lmultiflorum.tremona.fa)
REF_GFF=$(readlink -f reference_data/lmultiflorum.tremona.gene_annotation.liftoff.gff)
OUTDIR="results/annotation/"
T=4
ROOT=$(pwd)
run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

SUMMARY_FILE="mapping_summary.txt"
N_REF=$(awk -F'\t' '$3=="gene"' "${REF_GFF}" | wc -l)

for TARGET_FA in "${ROOT}"/reference_data/*.fa "${ROOT}"/reference_data/*.fasta; do
    [ -e "$TARGET_FA" ] || continue
    if [[ "$TARGET_FA" == *.fai ]]; then continue; fi
    if [[ "$TARGET_FA" == "$REF_FA" ]]; then continue; fi

    BASENAME=$(basename "${TARGET_FA}")
    PREFIX="${BASENAME%.*}"
    OUT_GFF="tremona_to_${PREFIX}.gff"

    if [[ -s "${OUT_GFF}" || -s "${OUT_GFF}.gz" ]]; then
        echo "[$(date)] Skipping ${PREFIX}, output already exists."
        continue
    fi
    
    run liftoff \
        -g "${REF_GFF}" \
        -o "tremona_to_${PREFIX}.gff" \
        -u "tremona_to_${PREFIX}.unmapped.txt" \
        -dir "liftoff_${PREFIX}_tmp" \
        -p "${T}" \
        -copies \
        -polish \
        -mm2_options="-a --end-bonus 5 --eqx -N 50 -p 0.5 -I 4G" \
        "${TARGET_FA}" "${REF_FA}"

    if [ -f "tremona_to_${PREFIX}.gff_polished" ]; then
        mv "tremona_to_${PREFIX}.gff_polished" "tremona_to_${PREFIX}.gff"
    fi

    rm -rf "liftoff_${PREFIX}_tmp"
done

run pigz -p "${T}" *.gff