#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
SNPS_LD="${ROOT}/results/filtered_vcf/snps.ld.vcf.gz"
POP="${ROOT}/scripts/pop_gen/pop.tsv"
PAL_POP="${ROOT}/scripts/pop_gen/palette_pop.tsv"
PAL_ANC="${ROOT}/scripts/pop_gen/palete_ancestral.tsv"
OUT="${ROOT}/results/pca"
RSCRIPT="${ROOT}/scripts/pop_gen/goal1_pca_structure.R"

mkdir -p "${OUT}"
run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

# recompute=FALSE : recharge le projet sNMF existant au lieu de tout recalculer
run Rscript "${RSCRIPT}" "${SNPS_LD}" "${POP}" "${OUT}" "chr" 2 \
    "${PAL_POP}" "${PAL_ANC}" 6 FALSE