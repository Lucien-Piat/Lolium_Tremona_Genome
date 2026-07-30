#!/bin/bash
#SBATCH --job-name=ragtag
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --time=08:00:00
#SBATCH --output=logs/04d_ragtag_%x_%j.log

# Usage: sbatch scripts/04d_ragtag.sh <reference.fa.gz> <assembly.fa.gz> <label>

set -euo pipefail

# Argument check
[[ $# -eq 3 ]] || { echo "Usage: sbatch $0 <reference.fa.gz> <assembly.fa.gz> <label>" >&2; exit 1; }

REF_GZ=$(readlink -f "$1")
ASM_GZ=$(readlink -f "$2")
LABEL="$3"

SIF=$(readlink -f images/sif/assembly_tools.sif)
OUTDIR="results/06_ragtag/${LABEL}"
FINAL="${OUTDIR}/lmultiflorum.${LABEL}.scaffolded.fa"
T=${SLURM_CPUS_PER_TASK:-4}
BIND="/cluster/scratch"

mkdir -p "${OUTDIR}" logs

# File existence checks
for f in "${SIF}" "${ASM_GZ}" "${REF_GZ}"; do
    [[ -s "$f" ]] || { echo "ERROR: Missing $f" >&2; exit 1; }
done

run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }

[[ -s "${OUTDIR}/ref.fa" ]] || run pigz -dcp "${T}" "${REF_GZ}" > "${OUTDIR}/ref.fa"
[[ -s "${OUTDIR}/query.fa" ]] || run pigz -dcp "${T}" "${ASM_GZ}" > "${OUTDIR}/query.fa"

# Run RagTag
SCAFFOLD="${OUTDIR}/ragtag_out/ragtag.scaffold.fasta"
if [[ ! -s "${SCAFFOLD}" ]]; then
    run ragtag.py scaffold \
        "${OUTDIR}/ref.fa" \
        "${OUTDIR}/query.fa" \
        -o "${OUTDIR}/ragtag_out" \
        -t "${T}" -u -r --aligner minimap2
fi

[[ -s "${SCAFFOLD}" ]] || { echo "ERROR: RagTag failed to produce ${SCAFFOLD}" >&2; exit 1; }

# Finalize and compress
mv "${SCAFFOLD}" "${FINAL}"
run pigz -p "${T}" "${FINAL}"

# Cleanup intermediate fasta files
rm -f "${OUTDIR}/ref.fa" "${OUTDIR}/query.fa"