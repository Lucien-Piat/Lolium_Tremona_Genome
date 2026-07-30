#!/bin/bash
#SBATCH --job-name=g5_pixy
#SBATCH --array=1-7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512M
#SBATCH --time=15:00:00
#SBATCH --output=logs/g5_%x_%A_%a.log
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
ALLSITES="${ROOT}/results/09_filtered_vcf/allsites.hc.vcf.gz"
POP="${ROOT}/scripts/09_pop_gen/pop.tsv"
FAI="${ROOT}/reference_data/lmultiflorum.tremona.fa.fai"
OUT="${ROOT}/results/09_diversity"

WINDOW=20000

export APPTAINER_HOME="${ROOT}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}" "${OUT}"

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FAI}" | cut -f1)
echo "[$(date)] tache ${SLURM_ARRAY_TASK_ID} -> ${CHR}"

run pixy --stats pi watterson_theta tajima_d \
    --vcf "${ALLSITES}" \
    --populations "${POP}" \
    --window_size "${WINDOW}" \
    --chromosomes "${CHR}" \
    --n_cores 1 \
    --output_folder "${OUT}" \
    --output_prefix "pixy_${CHR}"

