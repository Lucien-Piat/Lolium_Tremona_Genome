#!/bin/bash
#SBATCH --job-name=longqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=3G
#SBATCH --time=02:00:00
#SBATCH --output=logs/01_longqc_%j.log

set -euo pipefail

SIF="images/sif/qc_tools.sif"
READS="raw_reads/lmultiflorum_hifi.fastq.gz"
OUTDIR="results/01_qc/longqc"
T=${SLURM_CPUS_PER_TASK}

mkdir -p "$(dirname "${OUTDIR}")" logs
rm -rf "${OUTDIR}"

singularity exec "${SIF}" \
    longQC.py sampleqc -x pb-hifi -p "${T}" -o "${OUTDIR}" "${READS}"