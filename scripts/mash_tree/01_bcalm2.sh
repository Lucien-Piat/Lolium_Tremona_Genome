#!/bin/bash
# Usage:
#   1) Generate sample list:  ls $FILT_DIR/*_1.fastq.gz | sed 's/.*\///;s/_1.fastq.gz//' > samples.txt
#   2) Submit:                sbatch --array=1-$(wc -l < samples.txt) 01_bcalm2.sh
#SBATCH --job-name=bcalm2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=3G
#SBATCH --time=01:30:00
#SBATCH --output=logs/bcalm2_%A_%a.log

set -euo pipefail
exec 2>&1

# Configurable paths
SIF="${SIF_MASHTREE:-./image/mashtree_fastp_bcalm.sif}"
FILT_DIR="${FILT_DIR:-./filtered_reads}"
OUTDIR="${ASSEMBLIES_DIR:-./assemblies}"
SAMPLE_LIST="${SAMPLE_LIST:-./samples.txt}"
THREADS=6

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
if [[ -z "$SAMPLE" ]]; then
    echo "No sample at index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

WORKDIR="./tmp_${SAMPLE}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

mkdir -p "$OUTDIR"

if [[ -f "${OUTDIR}/${SAMPLE}.unitigs.fa.gz" ]]; then
    log "[SKIP] $SAMPLE already assembled"
    exit 0
fi

if ! mkdir "$WORKDIR" 2>/dev/null; then
    log "[SKIP] $SAMPLE already being processed"
    exit 0
fi

RUN="apptainer exec --no-home --pwd $(pwd) --bind $(pwd):$(pwd) $SIF"

cleanup() { rm -rf "$WORKDIR"; }
trap cleanup EXIT

log "[BCALM2] $SAMPLE (task ${SLURM_ARRAY_TASK_ID})"

echo "$(pwd)/${FILT_DIR}/${SAMPLE}_1.fastq.gz" >  "${WORKDIR}/input.fof"
echo "$(pwd)/${FILT_DIR}/${SAMPLE}_2.fastq.gz" >> "${WORKDIR}/input.fof"

$RUN bcalm \
    -in "${WORKDIR}/input.fof" \
    -out "${WORKDIR}/${SAMPLE}" \
    -kmer-size 31 \
    -abundance-min 2 \
    -nb-cores "$THREADS"

$RUN pigz -p "$THREADS" "${WORKDIR}/${SAMPLE}.unitigs.fa"
mv "${WORKDIR}/${SAMPLE}.unitigs.fa.gz" "${OUTDIR}/"

log "[DONE] $SAMPLE -> ${OUTDIR}/${SAMPLE}.unitigs.fa.gz ($(du -h "${OUTDIR}/${SAMPLE}.unitigs.fa.gz" | cut -f1))"