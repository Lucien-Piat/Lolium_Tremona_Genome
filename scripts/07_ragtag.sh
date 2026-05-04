#!/bin/bash
#SBATCH --job-name=ragtag
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=05:00:00
#SBATCH --output=logs/04d_ragtag_%x_%j.log

# Usage: sbatch scripts/04d_ragtag.sh <assembly.fa.gz> <label>

set -euo pipefail

ASM_GZ=$(readlink -f "$1")
LABEL="$2"

SIF=$(readlink -f images/sif/evaluation.sif)
REF_GZ=$(readlink -f reference_data/ciao_unp.fasta.gz)
OUTDIR="results/04d_ragtag/${LABEL}"
T=${SLURM_CPUS_PER_TASK:-4}

BIND="/cluster/scratch"

run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }

mkdir -p "${OUTDIR}" logs

echo "Label:      ${LABEL}"
echo "Assembly:   ${ASM_GZ}"
echo "Reference:  ${REF_GZ}"
echo "Threads:    ${T}"
echo "Started:    $(date)"

run pigz -dcp "${T}" "${REF_GZ}" > "${OUTDIR}/ref.fa"
run pigz -dcp "${T}" "${ASM_GZ}" > "${OUTDIR}/query.fa"

run ragtag.py scaffold \
    "${OUTDIR}/ref.fa" \
    "${OUTDIR}/query.fa" \
    -o "${OUTDIR}/ragtag_out" \
    -t "${T}" \
    -u \
    -r \
    --aligner minimap2

SCAFFOLD="${OUTDIR}/ragtag_out/ragtag.scaffold.fasta"
[[ -f "${SCAFFOLD}" ]] || { echo "ERROR: no RagTag output" >&2; exit 0; }

FINAL="${OUTDIR}/lmultiflorum.${LABEL}.scaffolded.fa"
mv "${SCAFFOLD}" "${FINAL}"
run pigz -p "${T}" "${FINAL}"

rm -f "${OUTDIR}/ref.fa" "${OUTDIR}/query.fa"

echo "Done at $(date)"
echo "Output: ${FINAL}.gz"