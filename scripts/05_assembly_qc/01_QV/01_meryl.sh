#!/bin/bash
#SBATCH --job-name=meryl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=4:00:00
#SBATCH --output=logs/02_meryl_%j.log

# Boilerplate
set -euo pipefail
SIF="images/sif/read_profiling.sif"
READS="./lmultiflorum_tremona_hifi.fastq.gz"
MERYL_DIR="results/meryl"
T=${SLURM_CPUS_PER_TASK}
mkdir -p "${MERYL_DIR}" logs

# Build 21mer db
singularity exec "${SIF}" \
    meryl count k=21 threads=${T} memory=150 "${READS}" output "${MERYL_DIR}/lmultiflorum.meryl"

# Build the histogram
singularity exec "${SIF}" \
    meryl histogram "${MERYL_DIR}/lmultiflorum.meryl" > "${MERYL_DIR}/lmultiflorum.hist"

if tar cf - -C "${MERYL_DIR}" lmultiflorum.meryl \
     | singularity exec "${SIF}" pigz -p "${T}" \
     > "${MERYL_DIR}/lmultiflorum.meryl.tar.gz"; then
    rm -rf "${MERYL_DIR}/lmultiflorum.meryl"
else
    echo 'ERROR' >&2
    exit 1
fi