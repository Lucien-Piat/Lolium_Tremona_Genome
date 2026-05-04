#!/bin/bash
#SBATCH --job-name=pg_main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --time=24:00:00
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
mkdir -p "${CACHE}/matplotlib"

run() {
    singularity exec --bind "${BIND}" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

BUSCO_TABLE=$(find "${OUTDIR}/compleasm_out" -name 'full_table_busco_format.tsv' | head -1)
[[ -n "${BUSCO_TABLE}" && -f "${BUSCO_TABLE}" ]] || \
    { echo "ERROR: full_table_busco_format.tsv not found, did compleasm finish?" >&2; exit 1; }

ASM_FA="${OUTDIR}/assembly.fa"
ASM_FAI="${OUTDIR}/assembly.fa.fai"
PURGEHAP_LOG="${OUTDIR}/curated.contig_associations.log"

for f in "${ASM_FA}" "${ASM_FAI}" "${PURGEHAP_LOG}" "${BUSCO_TABLE}" "${TRANSCRIPTS}"; do
    [[ -f "${f}" ]] || { echo "ERROR: missing input: ${f}" >&2; exit 1; }
done

cd "${OUTDIR}"


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

FINAL_RAW="final_primary_with_trim.fa"
FINAL_OUT="lmultiflorum.purgegrass.fa"

if [[ -s "${FINAL_RAW}" ]]; then
    mv "${FINAL_RAW}" "${FINAL_OUT}"
    [[ -f "${FINAL_RAW}.stats" ]] && mv "${FINAL_RAW}.stats" "${FINAL_OUT}.stats"
    run pigz -p "${T}" "${FINAL_OUT}"
    echo "Final purged assembly: ${OUTDIR}/${FINAL_OUT}.gz"
else
    echo "ERROR: expected ${FINAL_RAW} not found" >&2
    echo "Contents of ${OUTDIR}:"
    ls -lh
    exit 1
fi

cd "${ROOT}"
[[ -f "${TRACKING}" ]] || printf 'stage\tfile\tcontigs\tsize\n' > "${TRACKING}"
stats=$(run seqkit stats -T "${OUTDIR}/${FINAL_OUT}.gz" | tail -1)
nseq=$(echo "${stats}" | cut -f4)
size=$(echo "${stats}" | cut -f5)
printf 'assembly-purgegrass\t%s\t%s\t%s\n' \
    "$(readlink -f "${OUTDIR}/${FINAL_OUT}.gz")" "${nseq}" "${size}" >> "${TRACKING}"
echo "assembly-purgegrass: ${nseq} contigs, ${size} bp"

rm -rf "${CACHE}"