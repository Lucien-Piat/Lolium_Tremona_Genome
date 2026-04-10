#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3G
#SBATCH --time=24:00:00
#SBATCH --output=logs/03_hifiasm_%j.log

SIF="images/sif/hifiasm.sif"
READS="raw_reads/lmultiflorum_hifi.fastq.gz"
OUTDIR="results/02_assembly"
PREFIX="${OUTDIR}/lmultiflorum"
THREADS=${SLURM_CPUS_PER_TASK}

mkdir -p "${OUTDIR}" logs

singularity exec "${SIF}" \
    hifiasm \
        -t "${THREADS}" \
        -o "${PREFIX}" \
        --telo-m TTTAGGG \
        "${READS}"

for TYPE in bp.p_ctg bp.hap1.p_ctg bp.hap2.p_ctg; do
    awk '/^S/{print ">"$2; print $3}' "${PREFIX}.${TYPE}.gfa" \
        | singularity exec "${SIF}" pigz -p "${THREADS}" \
        > "${PREFIX}.${TYPE}.fa.gz"
done

rm -f "${PREFIX}"*.bin

for GFA in "${PREFIX}"*.gfa; do
    singularity exec "${SIF}" pigz -p "${THREADS}" "${GFA}"
done