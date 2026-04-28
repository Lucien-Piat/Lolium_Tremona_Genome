#!/bin/bash
#SBATCH --job-name=redundans
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=18:00:00
#SBATCH --output=logs/04_redundans_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/redundans.sif)
ASM_GZ=$(readlink -f results/02_assembly/lmultiflorum.bp.p_ctg.fa.gz)
OUTDIR="results/04_redundans"
TRACKING="results/assembly_tracking.tsv"
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)

BIND="/cluster/scratch"

CACHE="$(pwd)/.cache_redundans"
mkdir -p "${CACHE}/home" "${CACHE}/matplotlib"

run() {
    singularity exec \
        --bind "${BIND}" \
        --env HOME="${CACHE}/home" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

track() {
    local stage=$1 fa=$2
    local stats nseq size
    stats=$(singularity exec --bind "${BIND}" images/sif/polish.sif seqkit stats -T "${fa}" | tail -1)
    nseq=$(echo "${stats}" | cut -f4)
    size=$(echo "${stats}" | cut -f5)
    printf '%s\t%s\t%s\t%s\n' "${stage}" "$(readlink -f "${fa}")" "${nseq}" "${size}" >> "${TRACKING}"
    echo "${stage}: ${nseq} contigs, ${size} bp"
}

mkdir -p "${OUTDIR}" logs
[[ -f "${TRACKING}" ]] || printf 'stage\tfile\tcontigs\tsize\n' > "${TRACKING}"

[[ -f "${ASM_GZ}" ]] || { echo "ERROR: input missing: ${ASM_GZ}" >&2; exit 1; }

ASM_FA="${OUTDIR}/input_assembly.fa"
if [[ ! -f "${ASM_FA}" ]]; then
    echo "Decompressing assembly..."
    singularity exec --bind "${BIND}" images/sif/polish.sif \
        pigz -dcp "${T}" "${ASM_GZ}" > "${ASM_FA}"
fi

RUNDIR="${OUTDIR}/run"
rm -rf "${RUNDIR}"

run redundans.py \
    -v \
    -f "${ASM_FA}" \
    -o "${RUNDIR}" \
    -t "${T}" \
    --noscaffolding \
    --nogapclosing


REDUCED="${RUNDIR}/contigs.reduced.fa"
[[ -s "${REDUCED}" ]] || { echo "ERROR: redundans did not produce reduced contigs at ${REDUCED}" >&2; exit 1; }

cp "${REDUCED}" "${OUTDIR}/lmultiflorum.reduced.fa"
singularity exec --bind "${BIND}" images/sif/polish.sif \
    pigz -p "${T}" "${OUTDIR}/lmultiflorum.reduced.fa"


rm -f "${ASM_FA}"
rm -rf "${CACHE}"


cd "${ROOT}"
track "assembly-reduced" "${OUTDIR}/lmultiflorum.reduced.fa.gz"
