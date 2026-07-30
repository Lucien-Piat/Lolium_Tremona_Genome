#!/bin/bash
#SBATCH --job-name=te_rm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=712M
#SBATCH --time=48:00:00
#SBATCH --output=logs/05_04_rm_%j.log

# Usage : sbatch scripts/07_te_annotation/04_repeatmasker.sh <chr>

set -euo pipefail
SIF=$(readlink -f images/sif/hite.sif)
CHUNKDIR=$(readlink -f results/07_te_hite/chunks)
LIB=$(readlink -f results/07_te_hite/library/tremona_TE.lib.fa)
OUTBASE="results/07_te_hite/annotation"
CHR=${1:?usage: sbatch 04_repeatmasker.sh <chr>}
T=${SLURM_CPUS_PER_TASK:-16}
PA=$(( T / 4 )); [[ "${PA}" -ge 1 ]] || PA=1
BIND=$PWD
: "${TMPDIR:=/tmp}"
run() { singularity exec --bind "${BIND}" --bind "${TMPDIR}" "${SIF}" "$@"; }
mkdir -p "${OUTBASE}" logs

CHUNK=$(readlink -f "${CHUNKDIR}/${CHR}.fa")
OUTDIR="${OUTBASE}/${CHR}"
mkdir -p "${OUTDIR}"
[[ -s "${CHUNK}" && -s "${LIB}" ]] || { echo "ERROR: chunk ou librairie manquant" >&2; exit 1; }

run RepeatMasker \
    -pa "${PA}" \
    -lib "${LIB}" \
    -a \
    -xsmall \
    -gff \
    -no_is \
    -dir "${OUTDIR}" \
    "${CHUNK}"

echo "${CHR} masque -> ${OUTDIR}"