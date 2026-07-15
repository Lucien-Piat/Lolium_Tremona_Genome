#!/bin/bash
#SBATCH --job-name=g6_dxy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=08:00:00
#SBATCH --output=logs/g6_%x_%j.log
set -euo pipefail
# Goal 6 : Dxy par paire d'accessions CH (chaque CH = une population),

ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
ALLSITES="${ROOT}/results/filtered_vcf/allsites.hc.vcf.gz"
POP_CH="${ROOT}/scripts/pop_gen/pop_ch.tsv"
OUT="${ROOT}/results/dxy_ch"
WINDOW=20000

export APPTAINER_HOME="${ROOT}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}" "${OUT}"

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

if [ ! -s "${OUT}/ch_only.vcf.gz" ]; then
    run bash -c "cut -f1 '${POP_CH}' > '${OUT}/ch_samples.txt'"
    run bash -c "bcftools view -S '${OUT}/ch_samples.txt' '${ALLSITES}' -Oz -o '${OUT}/ch_only.vcf.gz'"
    run tabix -p vcf "${OUT}/ch_only.vcf.gz"
fi

run pixy --stats dxy \
    --vcf "${OUT}/ch_only.vcf.gz" \
    --populations "${POP_CH}" \
    --window_size "${WINDOW}" \
    --n_cores 4 \
    --output_folder "${OUT}" \
    --output_prefix "ch"
