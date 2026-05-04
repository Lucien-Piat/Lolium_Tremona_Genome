#!/bin/bash
#SBATCH --job-name=pg_purgehap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --time=16:00:00
#SBATCH --output=logs/04c_pg_02_purgehap_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/purgegrass.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
OUTDIR="results/04c_purgegrass"
T=${SLURM_CPUS_PER_TASK:-4}

BIND="/cluster/scratch"

CACHE="$(pwd)/.cache_purgegrass"
mkdir -p "${CACHE}/home"

run() {
    singularity exec --bind "${BIND}" "${SIF}" "$@"
}
run_sh() {
    singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"
}

cd "${OUTDIR}"

if [[ ! -f "aligned.bam" ]]; then
    echo "Step 1: minimap2 mapping at $(date)"
    run_sh "minimap2 -t 12 -I 5G -K 5G -ax map-hifi --secondary=no \
              assembly.fa ${READS} \
            | samtools sort -m 2G -@ 4 -o aligned.bam -"
    run samtools index -@ 4 aligned.bam
else
    echo "Step 1: reusing existing aligned.bam"
fi

GENCOV=$(ls aligned.bam*gencov 2>/dev/null | head -1 || true)
if [[ -z "${GENCOV}" ]]; then
    echo "Step 2: purge_haplotigs hist at $(date)"
    run purge_haplotigs hist -b aligned.bam -g assembly.fa -t 8
    GENCOV=$(ls aligned.bam*gencov 2>/dev/null | head -1)
else
    echo "Step 2: reusing existing ${GENCOV}"
fi

if [[ -z "${GENCOV}" ]] || [[ ! -f "${GENCOV}" ]]; then
    echo "ERROR: no gencov file produced by purge_haplotigs hist" >&2
    ls -lh aligned.bam* >&2
    exit 1
fi
echo "Using gencov file: ${GENCOV}"

run purge_haplotigs cov -i "${GENCOV}" -l 5 -m 40 -h 100 -j 101 -s 80

run purge_haplotigs purge -t "${T}" -g assembly.fa -c coverage_stats.csv -a 70 -I 5G

rm -rf "${CACHE}"