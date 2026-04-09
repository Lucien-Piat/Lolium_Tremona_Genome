#!/bin/bash
#SBATCH --job-name=longqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=logs/01_longqc_%j.log

SIF="images/sif/qc_tools.sif"
READS="raw_reads/lmultiflorum_hifi.fastq.gz"
OUTDIR="results/01_qc/longqc"
THREADS=${SLURM_CPUS_PER_TASK}

mkdir -p "${OUTDIR}" logs

singularity exec "${SIF}" \
    longQC.py sampleqc \
        -x pb-hifi \
        -o "${OUTDIR}" \
        -t "${THREADS}" \
        "${READS}"