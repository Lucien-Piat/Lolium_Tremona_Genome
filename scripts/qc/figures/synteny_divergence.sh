#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
DATA="${ROOT}/results/data_circo"
run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

run python3 scripts/qc/circo_plot/lib/synteny_divergence.py \
    --collinearity "${DATA}/mcscanx/tremona.collinearity" \
    --blast        "${DATA}/mcscanx/tremona.blast" \
    --output       "${DATA}/synteny_divergence.pdf"