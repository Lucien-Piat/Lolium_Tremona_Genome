#!/bin/bash
#SBATCH --job-name=hite_dn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --output=logs/05_02_hite_%j.log

set -euo pipefail
SIF=$(readlink -f images/sif/hite.sif)
CHUNKDIR=$(readlink -f results/te_hite/chunks)
OUTBASE="results/te_hite/denovo"
CHR=${1:?usage: sbatch 02_hite_denovo.sh <chr_name>}
T=${SLURM_CPUS_PER_TASK:-32}
BIND=$PWD

: "${TMPDIR:=/tmp}"
mkdir -p "${TMPDIR}/mpl"
export SINGULARITYENV_TMPDIR="${TMPDIR}"
export SINGULARITYENV_TMP="${TMPDIR}"
export SINGULARITYENV_MPLCONFIGDIR="${TMPDIR}/mpl"
run() { singularity exec --bind "${BIND}" --bind "${TMPDIR}" "${SIF}" "$@"; }
mkdir -p "${OUTBASE}" logs

CHUNK=$(readlink -f "${CHUNKDIR}/${CHR}.fa")
OUTDIR=$(readlink -f "${OUTBASE}")/${CHR}
mkdir -p "${OUTDIR}"
[[ -s "${CHUNK}" ]] || { echo "ERROR: ${CHUNK} introuvable" >&2; exit 1; }

run python /HiTE/main.py \
    --genome "${CHUNK}" \
    --out_dir "${OUTDIR}" \
    --thread "${T}" \
    --plant 1 \
    --miu 1.3e-8 \
    --annotate 0 \
    --recover 1

LIB="${OUTDIR}/confident_TE.cons.fa"
[[ -s "${LIB}" ]] || { echo "ERROR: librairie vide pour ${CHR}" >&2; exit 1; }
echo "${CHR} done : $(grep -c '^>' "${LIB}") familles"