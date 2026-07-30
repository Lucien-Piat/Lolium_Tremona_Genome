#!/bin/bash
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/genome_analysis.sif"
LIB="${ROOT}/scripts/05_assembly_qc/02_circo_plot/lib"
DATA="${ROOT}/results/05_genome_landscape"

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

cd "${LIB}"

COMMON_ARGS=(
    --genome      "${ROOT}/reference_data/lmultiflorum.tremona.fa"
    --fai         "${ROOT}/reference_data/lmultiflorum.tremona.fa.fai"
    --gff         "${ROOT}/reference_data/lmultiflorum.tremona.gene_annotation.gff"
    --busco       "${ROOT}/reference_data/lmultiflorum.tremona_full_table_busco_format.tsv"
    --synteny     "${ROOT}/results/02_synteny/tremona_purged/self_synteny_links.tsv"
    --numt        "${DATA}/numt_links.tsv"
    --nupt        "${DATA}/nupt_links.tsv"
    --mito-fasta  "${DATA}/lmul_tremona.mito.fasta"
    --pltd-fasta  "${DATA}/lmul_tremona.pltd.fasta"
    --mito-gb     "${DATA}/lmul_tremona.mito.gb"
    --pltd-gb     "${DATA}/lmul_tremona.pltd.gb"
)

run python3 plot_composite.py "${COMMON_ARGS[@]}" \
    --te-gff      "${ROOT}/results/07_te_hite/tremona_TE.gff3" \
    --te-mapping  "${ROOT}/results/07_te_hite/tremona_TE.family_table.tsv" \
    --dist        "${ROOT}/results/07_te_hite/gene_te/te_gene_distance.tsv" \
    --classtab    "${ROOT}/results/07_te_hite/tremona_TE.class_table.tsv" \
    --partition   "${ROOT}/results/07_te_hite/genome_partition.tsv" \
    --output      "${DATA}/circos_composite.pdf"