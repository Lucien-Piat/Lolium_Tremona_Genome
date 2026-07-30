#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
SNPS_LD="${ROOT}/results/09_filtered_vcf/snps.ld.vcf.gz"
POP="${ROOT}/scripts/09_pop_gen/pop.tsv"
PAL_POP="${ROOT}/scripts/09_pop_gen/palette_pop.tsv"
OUT="${ROOT}/results/09_njtree"
RSCRIPT="${ROOT}/scripts/09_pop_gen/03_njtree.R"

mkdir -p "${OUT}"
run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

run Rscript "${RSCRIPT}" "${SNPS_LD}" "${POP}" "${OUT}" "chr" "${PAL_POP}"