#!/bin/bash
#SBATCH --job-name=oatk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --output=logs/06_oatk_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/oatk.sif)
POLISH_SIF=$(readlink -f images/sif/polish.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
MITO_DB=$(readlink -f reference_data/oatkdb/OatkDB/v20230921/embryophyta_mito.fam)
PLTD_DB=$(readlink -f reference_data/oatkdb/OatkDB/v20230921/embryophyta_pltd.fam)
OUTDIR="results/05_oatk"
PREFIX="lmultiflorum"
T=${SLURM_CPUS_PER_TASK}
ROOT=$(pwd)
BIND="/cluster/scratch"
run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_polish() { singularity exec --bind "${BIND}" "${POLISH_SIF}" "$@"; }
mkdir -p "${OUTDIR}" logs
cd "${OUTDIR}"

# Subsample 
run_polish bash -c "seqkit sample -s 42 -p 0.5 -j ${T} ${READS} -o reads_subsampled.fastq.gz"

SYNCMER_COV=150

# Organelle assembly
run oatk \
    -k 1001 \
    -c "${SYNCMER_COV}" \
    -t "${T}" \
    -m "${MITO_DB}" \
    -p "${PLTD_DB}" \
    -o "${PREFIX}" \
    reads_subsampled.fastq.gz