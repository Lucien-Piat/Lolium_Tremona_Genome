#!/bin/bash
# Run PurgeGrass.sh locally

set -uo pipefail

ROOT=$(pwd)
OUTDIR="${ROOT}/results/04c_purgegrass"
SIF="${ROOT}/images/sif/purgegrass.sif"
PURGEGRASS_DIR="${ROOT}/PhaseGrass/PurgeGrass"
TRANSCRIPTS="${ROOT}/reference_data/transcripts/lolmu_transcripts.fa"
T=8
BUSCO_TABLE="full_table_busco_format.tsv"

cd "${OUTDIR}"

singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" \
    bash "${PURGEGRASS_DIR}/PurgeGrass.sh" \
    -a "${OUTDIR}/assembly.fa" \
    -b "${OUTDIR}/${BUSCO_TABLE}" \
    -f "${OUTDIR}/assembly.fa.fai" \
    -g "${TRANSCRIPTS}" \
    -m /opt/conda/bin/MCScanX \
    -p "${OUTDIR}/curated.contig_associations.log" \
    -s "${PURGEGRASS_DIR}/scripts" \
    -t "${T}" 2>&1 | tee purgegrass_local.log

PG_EXIT=${PIPESTATUS[0]}