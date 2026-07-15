#!/bin/bash
#SBATCH --job-name=g8_fst
#SBATCH --array=1-7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512M
#SBATCH --time=10:00:00
#SBATCH --output=logs/g8_%x_%A_%a.log
set -euo pipefail
# FST par fenetre (pixy), un chromosome par tache.

ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
ALLSITES="${ROOT}/results/filtered_vcf/allsites.hc.vcf.gz"
POP="${ROOT}/scripts/pop_gen/pop.tsv"
FAI="${ROOT}/reference_data/lmultiflorum.tremona.fa.fai"
OUT="${ROOT}/results/fst_windows"
WINDOW=20000

export APPTAINER_HOME="${ROOT}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}" "${OUT}"

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FAI}" | cut -f1)
echo "[$(date)] ${CHR}"

run pixy --stats fst \
    --vcf "${ALLSITES}" \
    --populations "${POP}" \
    --window_size "${WINDOW}" \
    --chromosomes "${CHR}" \
    --fst_type wc \
    --n_cores 1 \
    --output_folder "${OUT}" \
    --output_prefix "pixy_${CHR}"

echo "[$(date)] ${CHR} termine"


#cd results/fst_windows
#head -1 pixy_chr1_fst.txt > pixy_fst.txt
#for c in chr1 chr2 chr3 chr4 chr5 chr6 chr7; do
#    tail -n +2 pixy_${c}_fst.txt >> pixy_fst.txt
#done