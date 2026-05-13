#!/bin/bash
#SBATCH --job-name=pg_purgehap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=90G
#SBATCH --time=12:00:00
#SBATCH --output=logs/04c_pg_02_purgehap_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/purgegrass.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
OUTDIR="results/04c_purgegrass"
T=${SLURM_CPUS_PER_TASK:-4}
BIND="/cluster/scratch"
CACHE="$(pwd)/.cache_purgegrass"
mkdir -p "${CACHE}/home"
run() { singularity exec --bind "${BIND}" "${SIF}" "$@";}
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@";}

cd "${OUTDIR}"

if [[ ! -f "aligned.bam" ]]; then
    run_sh "minimap2 -t $((T-2)) -I 5G -K 5G -ax map-hifi --secondary=no \
              assembly.fa ${READS} \
            | samtools sort -m 3G -@ 2 -o aligned.bam -"
    run samtools index -@ 4 aligned.bam
fi

GENCOV=$(ls aligned.bam*gencov 2>/dev/null | head -1 || true)
if [[ -z "${GENCOV}" ]]; then
    run purge_haplotigs hist -b aligned.bam -g assembly.fa -t "${T}"
    GENCOV=$(ls aligned.bam*gencov 2>/dev/null | head -1)
fi

if [[ ! -f coverage_stats.csv ]]; then
    run purge_haplotigs cov -i "${GENCOV}" -l 5 -m 40 -h 100 -j 101 -s 80
fi

run purge_haplotigs purge -t "${T}" -g assembly.fa -c coverage_stats.csv -a 70 -I 1G
