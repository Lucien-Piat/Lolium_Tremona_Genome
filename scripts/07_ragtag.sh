#!/bin/bash
#SBATCH --job-name=ragtag
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --time=04:00:00
#SBATCH --output=logs/07_ragtag_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/evaluation.sif)
ASM_GZ=$(readlink -f results/05_blobtoolkit/lmultiflorum.decontam.fa.gz)
REF_GZ=$(readlink -f reference_data/ciao.fasta.gz)
OUTDIR="results/06_scaffolding"
TRACKING="results/assembly_tracking.tsv"
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)

BIND="/cluster/scratch"

run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }

track() {
    local stage=$1 fa=$2
    local stats nseq size
    stats=$(run seqkit stats -T "${fa}" | tail -1)
    nseq=$(echo "${stats}" | cut -f4)
    size=$(echo "${stats}" | cut -f5)
    printf '%s\t%s\t%s\t%s\n' "${stage}" "$(readlink -f "${fa}")" "${nseq}" "${size}" >> "${TRACKING}"
    echo "${stage}: ${nseq} contigs, ${size} bp"
}

mkdir -p "${OUTDIR}" logs
[[ -f "${TRACKING}" ]] || printf 'stage\tfile\tcontigs\tsize\n' > "${TRACKING}"

run pigz -dcp "${T}" "${REF_GZ}" > "${OUTDIR}/ref.fa"
run pigz -dcp "${T}" "${ASM_GZ}" > "${OUTDIR}/query.fa"

echo "Running RagTag scaffold..."
run ragtag.py scaffold \
    "${OUTDIR}/ref.fa" \
    "${OUTDIR}/query.fa" \
    -o "${OUTDIR}/ragtag_ciao" \
    -t "${T}" \
    -u \
    -r \
    --aligner minimap2

SCAFFOLD="${OUTDIR}/ragtag_ciao/ragtag.scaffold.fasta"
[[ -f "${SCAFFOLD}" ]] || { echo "ERROR: no RagTag output" >&2; exit 1; }

mv "${SCAFFOLD}" "${OUTDIR}/lmultiflorum.scaffolded.fa"
run pigz -p "${T}" "${OUTDIR}/lmultiflorum.scaffolded.fa"
rm -f "${OUTDIR}/ref.fa" "${OUTDIR}/query.fa"

cd "${ROOT}"
track "assembly-scaffolded" "${OUTDIR}/lmultiflorum.scaffolded.fa.gz"

echo "Done. RagTag output at ${OUTDIR}/ragtag_ciao/"