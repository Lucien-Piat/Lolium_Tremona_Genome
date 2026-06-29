#!/bin/bash
#SBATCH --job-name=te_landscape
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/05_06_landscape_%j.log

# Usage : sbatch scripts/te_annotation/06_landscape.sh

set -euo pipefail
SIF=$(readlink -f images/sif/hite.sif)
ANNDIR=$(readlink -f results/te_hite/annotation)
GENOME=$(readlink -f reference_data/lmultiflorum.tremona.fa)
OUTDIR="results/te_hite/landscape"
BIND=$PWD
: "${TMPDIR:=/tmp}"

RMLIB=$(singularity exec --bind "${BIND}" "${SIF}" bash -c \
    'dirname "$(find /opt/conda -name SearchResult.pm 2>/dev/null | head -1)"')
[[ -n "${RMLIB}" && "${RMLIB}" != "." ]] \
    || { echo "ERROR: SearchResult.pm introuvable dans le conteneur" >&2; exit 1; }
export SINGULARITYENV_PERL5LIB="${RMLIB}"
echo "PERL5LIB (conteneur) -> ${RMLIB}"

run() { singularity exec --bind "${BIND}" --bind "${TMPDIR}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}" logs

ALIGN="${OUTDIR}/tremona.combined.align"
DIVSUM="${OUTDIR}/tremona.divsum"
HTML="${OUTDIR}/tremona_landscape.html"

shopt -s nullglob
files=( "${ANNDIR}"/chr*/*.align )
[[ ${#files[@]} -gt 0 ]] || { echo "ERROR: aucun .align (04 lance sans -a ?)" >&2; exit 1; }
cat "${files[@]}" > "${ALIGN}"
echo "Combine ${#files[@]} fichiers .align"

GSIZE=$(awk '{s+=$2} END{print s}' "${GENOME}.fai")

run calcDivergenceFromAlign.pl -s "${DIVSUM}" "${ALIGN}"
[[ -s "${DIVSUM}" ]] || { echo "ERROR: .divsum vide" >&2; exit 1; }

run bash -c "createRepeatLandscape.pl -div '${DIVSUM}' -g ${GSIZE} > '${HTML}'"

echo "Landscape : ${HTML}"
echo "Table brute (replot avec axe temps via T = K/(2mu)) : ${DIVSUM}"