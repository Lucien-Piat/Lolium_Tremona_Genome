#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=24:00:00
#SBATCH --output=logs/03_hifiasm_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/hifiasm.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
OUTDIR="results/03_assembly"
T=${SLURM_CPUS_PER_TASK}
BIND="/cluster/scratch"
run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}" logs
cd "${OUTDIR}"

# Assembly
run hifiasm -t "${T}" -o lmultiflorum --telo-m TTTAGGG --hg-size 2.2g --purge-max 80 -s 0.5 -O 3 "${READS}"

# Extract the contigs
for TYPE in bp.p_ctg bp.hap1.p_ctg bp.hap2.p_ctg; do
    awk '/^S/{print ">"$2; print $3}' "lmultiflorum.${TYPE}.gfa" \
        | run pigz -p "${T}" > "lmultiflorum.${TYPE}.fa.gz"
done

# Compress and cleanup
rm -f lmultiflorum*.bin
run pigz -p "${T}" lmultiflorum*.gfa
