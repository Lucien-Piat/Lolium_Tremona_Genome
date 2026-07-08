#!/bin/bash
#SBATCH --job-name=filt_as
#SBATCH --array=1-7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=18:00:00
#SBATCH --output=logs/filter_%x_%A_%a.log
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
ALLSITES="${ROOT}/results/04_joint_calling/cohort_allsites.vcf.gz"
MASK="${ROOT}/results/mask/accessible.bed"
FAI="${ROOT}/reference_data/lmultiflorum.tremona.fa.fai"
TMP="${ROOT}/results/filtered_vcf/tmp_allsites"

MIN_DP=10
MIN_GQ=30

export APPTAINER_HOME="${ROOT}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}" "${TMP}" logs

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FAI}" | cut -f1)
OUT_CHR="${TMP}/${CHR}.hc.vcf.gz"
BED_CHR="${TMP}/${CHR}.bed"

if [ -s "${OUT_CHR}" ]; then
    exit 0
fi

run bash -c "awk -v c='${CHR}' '\$1==c' '${MASK}' > '${BED_CHR}'"
if [ ! -s "${BED_CHR}" ]; then
    echo "aucune region accessible sur ${CHR}" >&2
    exit 0
fi

if run bash -c "bcftools view -h '${ALLSITES}' | grep -q 'FORMAT=<ID=GQ,'"; then
    GT_FILT="FMT/DP<=${MIN_DP} | FMT/GQ<=${MIN_GQ}"
else
    GT_FILT="FMT/DP<=${MIN_DP}"
fi

run bash -c "set -o pipefail; \
    bcftools view -R '${BED_CHR}' '${ALLSITES}' -Ou \
    | bcftools +setGT - -Ou -- -t q -i '${GT_FILT}' -n . \
    | bcftools view -Oz -o '${OUT_CHR}'"
run tabix -p vcf "${OUT_CHR}"
