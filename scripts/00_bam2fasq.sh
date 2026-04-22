#!/bin/bash
#SBATCH --job-name=bam2fq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --time=02:00:00
#SBATCH --output=logs/00_bam2fq_%j.log

set -euo pipefail

SIF="images/sif/qc_tools.sif"
BAM="TREM1.hifi_reads.bam"
OUT="raw_reads/lmultiflorum_hifi.fastq.gz"
T=${SLURM_CPUS_PER_TASK}

mkdir -p raw_reads logs

# Half threads to samtools decode, half to pigz compress.
# samtools fastq is I/O bound here, pigz is CPU bound.
singularity exec "${SIF}" bash -c "
    samtools fastq -@ $((T/2)) '${BAM}' \
        | pigz -p $((T/2)) > '${OUT}'
"

echo 'Read stats:'
singularity exec "${SIF}" seqkit stats -a -T "${OUT}"