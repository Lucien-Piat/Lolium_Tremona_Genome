#!/bin/bash
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
LIB="${ROOT}/scripts/qc/circo_plot/lib"
DATA="${ROOT}/results/data_circo"

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

cd "${LIB}"

COMMON_ARGS=(
    --genome      "${DATA}/lmultiflorum.tremona.placed.fa"
    --fai         "${DATA}/lmultiflorum.tremona.placed.fa.fai"
    --gff         "${DATA}/tremona.gene_annotation.placed.gff"
    --busco       "${DATA}/full_table_busco_format.tsv"
    --synteny     "${DATA}/mcscanx/self_synteny_links.tsv"
    --numt        "${DATA}/numt_links.tsv"
    --nupt        "${DATA}/nupt_links.tsv"
    --mito-fasta  "${DATA}/lmul_tremona.mito.fasta"
    --pltd-fasta  "${DATA}/lmul_tremona.pltd.fasta"
    --mito-gb     "${DATA}/lmul_tremona.mito.gb"
    --pltd-gb     "${DATA}/lmul_tremona.pltd.gb"
    --coverage    "${DATA}/"TREM*_circos_50000bp.txt
)

echo "[$(date)] Building FULL variant"
run python3 plot_circos.py "${COMMON_ARGS[@]}" \
    --variant full \
    --output "${DATA}/circos_full.pdf"

echo "[$(date)] Building CLEAN variant"
run python3 plot_circos.py "${COMMON_ARGS[@]}" \
    --variant clean \
    --output "${DATA}/circos_clean.pdf"