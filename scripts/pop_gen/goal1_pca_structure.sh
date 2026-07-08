#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
SNPS_LD="${ROOT}/results/filtered_vcf/snps.ld.vcf.gz"
POP="${ROOT}/scripts/pop_gen/pop.tsv"
OUT="${ROOT}/results/pca"
RSCRIPT="${ROOT}/scripts/pop_gen/goal1_pca_structure.R"

mkdir -p "${OUT}"
run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

run Rscript "${RSCRIPT}" "${SNPS_LD}" "${POP}" "${OUT}" "chr" 2