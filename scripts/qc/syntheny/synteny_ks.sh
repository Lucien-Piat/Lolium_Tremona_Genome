#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
GENOME="${ROOT}/results/data_circo/lmultiflorum.tremona.placed.fa"
GFF="${ROOT}/results/data_circo/tremona.gene_annotation.placed.gff"
DATA="${ROOT}/results/data_circo"
MCS="${DATA}/mcscanx"
LIB="${ROOT}/scripts/qc/circo_plot/lib"

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

if [ ! -s "${MCS}/cds.fa" ]; then
    echo "[$(date)] Extracting CDS"
    run gffread -x "${MCS}/cds.fa" -g "${GENOME}" "${GFF}"
fi

run python3 "${LIB}/synteny_ks.py" \
    --collinearity "${MCS}/tremona.collinearity" \
    --cds          "${MCS}/cds.fa" \
    --output       "${DATA}/synteny_ks.pdf" \
    --keep-tmp