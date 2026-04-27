#!/bin/bash
#SBATCH --job-name=meryl_gs2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=15G
#SBATCH --time=4:00:00
#SBATCH --output=logs/02_meryl_genomescope_%j.log

set -euo pipefail

SIF="images/sif/qc_tools.sif"
READS="raw_reads/lmultiflorum_hifi.fastq.gz"
MERYL_DIR="results/01_qc/meryl"
GS_DIR="results/01_qc/genomescope2"
T=${SLURM_CPUS_PER_TASK}
KMER=21
MERYL_MEM=150 

mkdir -p "${MERYL_DIR}" "${GS_DIR}" logs

singularity exec "${SIF}" \
    meryl count k=${KMER} threads=${T} memory=${MERYL_MEM} \
    "${READS}" output "${MERYL_DIR}/lmultiflorum.meryl"

singularity exec "${SIF}" \
    meryl histogram "${MERYL_DIR}/lmultiflorum.meryl" \
    > "${MERYL_DIR}/lmultiflorum.hist"

singularity exec "${SIF}" \
    genomescope2 -i "${MERYL_DIR}/lmultiflorum.hist" \
    -o "${GS_DIR}" -k ${KMER} -p 2 --name_prefix lmultiflorum

if tar cf - -C "${MERYL_DIR}" lmultiflorum.meryl \
     | singularity exec "${SIF}" pigz -p "${T}" \
     > "${MERYL_DIR}/lmultiflorum.meryl.tar.gz"; then
    rm -rf "${MERYL_DIR}/lmultiflorum.meryl"
else
    echo 'ERROR' >&2
    exit 1
fi