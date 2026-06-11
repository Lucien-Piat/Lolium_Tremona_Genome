#!/bin/bash
# Usage: sbatch 01_fastp.sh [raw_reads_dir]
#SBATCH --job-name=fastp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=2G
#SBATCH --time=3:00:00
#SBATCH --output=logs/fastp_%j.log

set -euo pipefail
exec 2>&1

RAW_DIR="${RAW_DIR:-${1:-./raw_reads}}"
FILT_DIR="${FILT_DIR:-./filtered_reads}"
QC_DIR="${QC_DIR:-./fastp_reports}"
SIF="${SIF_MASHTREE:-./image/mashtree_fastp_bcalm.sif}"
THREADS=6

mkdir -p "$FILT_DIR" "$QC_DIR" logs

RUN="apptainer exec --no-home --pwd $(pwd) --bind $(pwd):$(pwd) $SIF"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

PROCESSED=0
SKIPPED=0

for r1 in "${RAW_DIR}"/*_1.fastq.gz; do
    [[ -e "$r1" ]] || continue
    SAMPLE=$(basename "$r1" _1.fastq.gz)
    r2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

    [[ -f "$r2" ]] || { log "[WARN] No R2 for $SAMPLE, skipping"; continue; }
    [[ -f "${FILT_DIR}/${SAMPLE}_1.fastq.gz" ]] && { SKIPPED=$((SKIPPED + 1)); continue; }

    log "[FASTP] $SAMPLE ($((PROCESSED + SKIPPED + 1)))"

    $RUN fastp \
        --in1 "$r1" \
        --in2 "$r2" \
        --out1 "${FILT_DIR}/${SAMPLE}_1.fastq.gz" \
        --out2 "${FILT_DIR}/${SAMPLE}_2.fastq.gz" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread "$THREADS" \
        --json "${QC_DIR}/${SAMPLE}.fastp.json" \
        --html /dev/null

    PROCESSED=$((PROCESSED + 1))
    log "[OK] $SAMPLE"
done
