#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=10:00:00
#SBATCH --output=logs/merqury_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/genome_analysis.sif)
ASM=$(readlink -f reference_data/lmultiflorum.tremona.fa)
MERYL_TAR=$(readlink -f results/meryl/lmultiflorum.meryl.tar.gz)
OUTDIR="results/05_merqury"
T=${SLURM_CPUS_PER_TASK:-8}
BIND="/cluster/scratch"
run()    { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"; }

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

run_sh "pigz -dc -p ${T} '${MERYL_TAR}' | tar -xf -"
MERYL_DIR=$(ls -d *.meryl 2>/dev/null | head -1 || true)

ln -sf "${ASM}" .
run_sh "export OMP_NUM_THREADS=${T} && \
    merqury.sh ${MERYL_DIR} $(basename ${ASM}) tremona"

rm -rf "${MERYL_DIR}"

