#!/bin/bash
# Run PurgeGrass.sh locally

set -uo pipefail

ROOT=$(pwd)
OUTDIR="${ROOT}/results/04c_purgegrass"
SIF="${ROOT}/images/sif/purgegrass.sif"
PURGEGRASS_DIR="${ROOT}/PhaseGrass/PurgeGrass"
TRANSCRIPTS="${ROOT}/reference_data/transcripts/lolmu_transcripts.fa"
T=4
BUSCO_TABLE="${ROOT}/results/04c_purgegrass/busco_table.tsv"
cd "${OUTDIR}"

singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" \
    bash "${PURGEGRASS_DIR}/PurgeGrass.sh" \
    -a assembly.fa \
    -b "${BUSCO_TABLE}" \
    -f assembly.fa.fai \
    -g "${TRANSCRIPTS}" \
    -m /opt/conda/bin/MCScanX \ \
    -p curated.contig_associations.log \
    -s "${PURGEGRASS_DIR}/scripts" \
    -t "${T}" 2>&1 | tee purgegrass_local.log

PG_EXIT=${PIPESTATUS[0]}
