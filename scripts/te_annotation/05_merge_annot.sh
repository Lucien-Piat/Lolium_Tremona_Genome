#!/bin/bash
#SBATCH --job-name=te_annomerge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/05_05_annomerge_%j.log

# Usage : sbatch scripts/te_annotation/05_merge_annotation.sh

set -euo pipefail
SIF=$(readlink -f images/sif/hite.sif)
ANNDIR=$(readlink -f results/te_hite/annotation)
GENOME=$(readlink -f reference_data/lmultiflorum.tremona.fa)
OUTDIR="${ANNDIR}/merged"
BIND=$PWD
: "${TMPDIR:=/tmp}"
run() { singularity exec --bind "${BIND}" --bind "${TMPDIR}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}" logs

GFF="${OUTDIR}/tremona_TE.gff3"
OUT="${OUTDIR}/tremona_TE.out"
SOFT="${OUTDIR}/tremona.softmasked.fa"
shopt -s nullglob

# concat des .masked, dans l'ordre du .fai
: > "${SOFT}"
while read -r chr _; do
    m="${ANNDIR}/${chr}/${chr}.fa.masked"
    [[ -s "${m}" ]] || { echo "WARN: ${m} absent, chromosome saute" >&2; continue; }
    cat "${m}" >> "${SOFT}"
done < "${GENOME}.fai"
run samtools faidx "${SOFT}"

# GFF3 fusionne
cat "${ANNDIR}"/chr*/*.out.gff | grep -v '^#' | sort -k1,1 -k4,4n > "${OUTDIR}/.body.gff3"
{ echo "##gff-version 3"; cat "${OUTDIR}/.body.gff3"; } > "${GFF}"
rm -f "${OUTDIR}/.body.gff3"

# .out fusionne 
first=$(ls "${ANNDIR}"/chr*/*.out | head -1)
head -3 "${first}" > "${OUT}"
for f in "${ANNDIR}"/chr*/*.out; do tail -n +4 "${f}"; done >> "${OUT}"

awk '!/^>/ { t += length($0); s = $0; m += gsub(/[acgt]/, "", s) }
     END { printf "Genome masque : %.2f%% (%d / %d bp)\n", 100*m/t, m, t }' "${SOFT}"

echo "GFF3        : ${GFF}"
echo "Soft-masked : ${SOFT}"
