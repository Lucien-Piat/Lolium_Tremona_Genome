#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/merqury_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/QC.sif)
ASM=$(readlink -f reference_data/lmultiflorum.tremona.primary.fa)
MERYL_TAR=$(readlink -f results/02_qc/meryl/lmultiflorum.meryl.tar.gz)
OUTDIR="results/02_qc/merqury"
T=${SLURM_CPUS_PER_TASK:-8}
BIND="/cluster/scratch"
run()    { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"; }

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

tar -xzf "${MERYL_TAR}"
MERYL_DIR=$(ls -d *.meryl 2>/dev/null | head -1 || true)

ln -sf "${ASM}" .
run_sh "export OMP_NUM_THREADS=${T} && \
    merqury.sh ${MERYL_DIR} $(basename ${ASM}) tremona"

tar -czf "${MERYL_TAR}.new" "${MERYL_DIR}"
mv "${MERYL_TAR}.new" "${MERYL_TAR}"
rm -rf "${MERYL_DIR}"

