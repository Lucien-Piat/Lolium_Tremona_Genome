#!/bin/bash
#SBATCH --job-name=gm_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=18G
#SBATCH --time=06:00:00
#SBATCH --output=logs/mask_%x_%j.log
set -euo pipefail

SIF="$(readlink -f images/sif/popgen.sif)"
BIND="$PWD"
GENOME="reference_data/lmultiflorum.tremona.fa"
IDX="results/mask/genmap_index"

export APPTAINER_HOME="${PWD}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}" results/mask logs

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${BIND}" "${SIF}" "$@"; }

GB=$(awk '{s+=$2} END{printf "%.2f", s/1e9}' "${GENOME}.fai")

rm -rf "${IDX}"
run genmap index -F "${GENOME}" -I "${IDX}"