#!/bin/bash
#SBATCH --job-name=pigz_comp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=156M
#SBATCH --time=02:00:00
#SBATCH --output=logs/compress_%x_%j.log

set -euo pipefail

# Syntax corrected for variable expansion
FILE="${1:?usage: sbatch scripts/compress.sh <file.fastq>}"
T="${SLURM_CPUS_PER_TASK:-4}"
SIF="$(readlink -f images/sif/varcall.sif)"
BIND="$PWD"

echo "Compressing ${FILE} using ${T} threads..."
singularity exec --bind "${BIND}" "${SIF}" pigz --best --force -p "${T}" "${FILE}"
echo "Compression complete. Generated ${FILE}.gz"