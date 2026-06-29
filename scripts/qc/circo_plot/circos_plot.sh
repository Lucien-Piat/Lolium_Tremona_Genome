#!/bin/bash
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
LIB="${ROOT}/scripts/qc/circo_plot/lib"
DATA="${ROOT}/results/data_circo"

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

cd "${LIB}"

COMMON_ARGS=(
    --genome      "${ROOT}/reference_data/lmultiflorum.tremona.fa"
    --fai         "${ROOT}/reference_data/lmultiflorum.tremona.fa.fai"
    --gff         "${ROOT}/reference_data/lmultiflorum.tremona.gene_annotation.gff"
    --busco       "${ROOT}/reference_data/lmultiflorum.tremona_full_table_busco_format.tsv"
    --synteny     "${ROOT}/results/synteny/tremona_purged/self_synteny_links.tsv"
    --numt        "${DATA}/numt_links.tsv"
    --nupt        "${DATA}/nupt_links.tsv"
    --mito-fasta  "${DATA}/lmul_tremona.mito.fasta"
    --pltd-fasta  "${DATA}/lmul_tremona.pltd.fasta"
    --mito-gb     "${DATA}/lmul_tremona.mito.gb"
    --pltd-gb     "${DATA}/lmul_tremona.pltd.gb"
)

echo "[$(date)] Building Circos plot with TE and SNP tracks"
run python3 plot_circos.py "${COMMON_ARGS[@]}" \
    --te-gff      "${ROOT}/results/te_hite/tremona_TE.gff3" \
    --te-mapping  "${ROOT}/results/te_hite/tremona_TE.family_table.tsv" \
    --vcf         "NA" \
    --output      "${DATA}/circos_custom.pdf"