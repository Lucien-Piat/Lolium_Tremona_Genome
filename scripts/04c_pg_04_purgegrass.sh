#!/bin/bash
#SBATCH --job-name=pg_main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/04c_pg_04_purgegrass_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/purgegrass.sif)
OUTDIR="results/04c_purgegrass"
TRACKING="results/assembly_tracking.tsv"
TRANSCRIPTS=$(readlink -f reference_data/transcripts/lolmu_transcripts.fa)
PURGEGRASS_DIR=$(readlink -f PhaseGrass/PurgeGrass)
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)

BIND="/cluster/scratch"

CACHE="$(pwd)/.cache_pg_main"
mkdir -p "${CACHE}/home" "${CACHE}/matplotlib"

run() {
    singularity exec --bind "${BIND}" \
        --env HOME="${CACHE}/home" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

# Find the BUSCO full_table.tsv from step 3
BUSCO_TABLE=$(find "${OUTDIR}/busco_out" -name 'full_table.tsv' | head -1)
[[ -n "${BUSCO_TABLE}" && -f "${BUSCO_TABLE}" ]] || \
    { echo "ERROR: BUSCO full_table.tsv not found, did step 3 finish?" >&2; exit 1; }

# Sanity check all required inputs
ASM_FA="${OUTDIR}/assembly.fa"
ASM_FAI="${OUTDIR}/assembly.fa.fai"
PURGEHAP_LOG="${OUTDIR}/curated.contig_associations.log"

for f in "${ASM_FA}" "${ASM_FAI}" "${PURGEHAP_LOG}" "${BUSCO_TABLE}" "${TRANSCRIPTS}"; do
    [[ -f "${f}" ]] || { echo "ERROR: missing input: ${f}" >&2; exit 1; }
done

cd "${OUTDIR}"

echo "Running PurgeGrass.sh at $(date)"
echo "Assembly:        ${ASM_FA}"
echo "BUSCO table:     ${BUSCO_TABLE}"
echo "FAI index:       ${ASM_FAI}"
echo "Transcripts:     ${TRANSCRIPTS}"
echo "MCScanX:         /opt/MCScanX"
echo "PurgeHap log:    ${PURGEHAP_LOG}"
echo "Scripts dir:     ${PURGEGRASS_DIR}/scripts"
echo "Threads:         ${T}"

run bash "${PURGEGRASS_DIR}/PurgeGrass.sh" \
    -a assembly.fa \
    -b "${BUSCO_TABLE}" \
    -f assembly.fa.fai \
    -g "${TRANSCRIPTS}" \
    -m /opt/MCScanX \
    -p curated.contig_associations.log \
    -s "${PURGEGRASS_DIR}/scripts" \
    -t "${T}"

echo "PurgeGrass.sh done at $(date)"

# The wrapper writes final_primary_with_trim.fa as its end product.
# Rename to project convention and compress.
FINAL_RAW="final_primary_with_trim.fa"
FINAL_OUT="lmultiflorum.purgegrass.fa"

if [[ -s "${FINAL_RAW}" ]]; then
    mv "${FINAL_RAW}" "${FINAL_OUT}"
    [[ -f "${FI