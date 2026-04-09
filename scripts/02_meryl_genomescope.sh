#!/bin/bash
#SBATCH --job-name=meryl_gs2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --time=06:00:00
#SBATCH --output=logs/02_meryl_genomescope_%j.log

SIF="images/sif/qc_tools.sif"
READS="raw_reads/lmultiflorum_hifi.fastq.gz"
MERYL_DIR="results/01_qc/meryl"
GS_DIR="results/01_qc/genomescope2"
THREADS=${SLURM_CPUS_PER_TASK}
KMER=21

mkdir -p "${MERYL_DIR}" "${GS_DIR}" logs

singularity exec "${SIF}" \
    meryl count \
        k=${KMER} \
        threads=${THREADS} \
        memory=64 \
        "${READS}" \
        output "${MERYL_DIR}/lmultiflorum.meryl"

singularity exec "${SIF}" \
    meryl histogram \
        "${MERYL_DIR}/lmultiflorum.meryl" \
        > "${MERYL_DIR}/lmultiflorum.hist"

# GenomeScope2 (diploid, p=2)
singularity exec "${SIF}" \
    genomescope2 \
        -i "${MERYL_DIR}/lmultiflorum.hist" \
        -o "${GS_DIR}" \
        -k ${KMER} \
        -p 2 \
        --name_prefix lmultiflorum

# Tar the meryl DB (cluster TOS)
tar cf - -C "${MERYL_DIR}" lmultiflorum.meryl \
    | singularity exec "${SIF}" pigz -p "${THREADS}" \
    > "${MERYL_DIR}/lmultiflorum.meryl.tar.gz"
rm -rf "${MERYL_DIR}/lmultiflorum.meryl"