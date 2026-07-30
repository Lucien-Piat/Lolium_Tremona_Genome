#!/bin/bash
#SBATCH --job-name=organellar_depth
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --output=logs/07_organellar_depth_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/assembly_tools.sif)
ASM=$(readlink -f results/04c_purgegrass/final_primary_with_trim.fa)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
IDS=$(readlink -f results/04c_purgegrass/p_organelar.txt)
OUTDIR="results/07_organellar_depth"
T=${SLURM_CPUS_PER_TASK}
BIND="/cluster/scratch"
run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}" logs
cd "${OUTDIR}"

# Extract the flagged contig
run seqkit grep -f "${IDS}" "${ASM}" > organellar_subset.fa
run samtools faidx organellar_subset.fa

# Map HiFi reads to the subset, sort, index
run bash -c "minimap2 -ax map-hifi -t ${T} organellar_subset.fa ${READS} \
    | samtools sort -@ ${T} -o aligned.bam -"
run samtools index aligned.bam

run samtools coverage aligned.bam > per_contig_coverage.tsv