#!/bin/bash
#SBATCH --job-name=te_merge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/05_03_merge_%j.log

set -euo pipefail
SIF=$(readlink -f images/sif/hite.sif)
DENOVO=$(readlink -f results/07_te_hite/denovo)
OUTDIR="results/07_te_hite/library"
T=${SLURM_CPUS_PER_TASK:-8}
BIND=$PWD
run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}" logs

RAW="${OUTDIR}/all_chr.raw.fa"
NR="${OUTDIR}/tremona_TE.lib.fa"

: > "${RAW}"
for d in "${DENOVO}"/*/; do
    c=$(basename "${d}")
    lib="${d%/}/confident_TE.cons.fa"
    [[ -s "${lib}" ]] || { echo "WARN: pas de librairie dans ${d}" >&2; continue; }
    sed "s/^>/>${c}__/" "${lib}" >> "${RAW}"
done

run cd-hit-est -i "${RAW}" -o "${NR}" \
    -c 0.8 -n 5 -aS 0.95 -aL 0.95 -M 0 -T "${T}" -d 0 -g 1
