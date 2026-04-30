#!/bin/bash
#SBATCH --job-name=bam2fq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=384M
#SBATCH --time=03:00:00
#SBATCH --output=logs/00_bam2fq_%j.log

set -euo pipefail

SIF="images/sif/qc_tools.sif"
BAM="TREM1.hifi_reads.bam"
OUT="raw_reads/lmultiflorum_hifi.fastq.gz"
T=${SLURM_CPUS_PER_TASK}

mkdir -p raw_reads logs

singularity exec "${SIF}" bash -c "samtools fastq -@ 1 '${BAM}' | pigz -p $((T-1)) > '${OUT}'"
