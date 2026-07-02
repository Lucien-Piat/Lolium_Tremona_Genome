#!/bin/bash
#SBATCH --job-name=gm_map
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=10:00:00
#SBATCH --output=logs/mask_%x_%j.log
set -euo pipefail

SIF="$(readlink -f images/sif/popgen.sif)"
BIND="$PWD"
T="${SLURM_CPUS_PER_TASK:-16}"
IDX="results/mask/genmap_index"
OUT="results/mask"

READLEN=150
ERRORS=2

export APPTAINER_HOME="${PWD}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}"

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${BIND}" "${SIF}" "$@"; }

echo "[$(date)] map K=${READLEN} E=${ERRORS} T=${T}"
run genmap map -K "${READLEN}" -E "${ERRORS}" \
    -I "${IDX}" -O "${OUT}/genmap" -bg -T "${T}"

run bash -c "awk '\$4>=1' '${OUT}/genmap.bedgraph' | bedtools merge -i - > '${OUT}/mask_mappable.bed'"
echo "[$(date)] mappable -> ${OUT}/mask_mappable.bed"