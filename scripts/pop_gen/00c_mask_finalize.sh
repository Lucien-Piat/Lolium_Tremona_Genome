#!/bin/bash
#SBATCH --job-name=gm_final
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512M
#SBATCH --time=01:00:00
#SBATCH --output=logs/mask_%x_%j.log
set -euo pipefail
# Masque d'accessibilite a DEUX couches (aligne sur Table 13) :
#   mappabilite (GenMap) INTER (non-TE).
# La couche profondeur est ecartee : les collapses de paralogues ont deja
# ete retires physiquement lors du self-synteny purge de l'assemblage.

SIF="$(readlink -f images/sif/popgen.sif)"
BIND="$PWD"
GENOME="reference_data/lmultiflorum.tremona.fa"
TE_GFF="results/te_hite/tremona_TE.gff3"
OUT="results/mask"
TMP="${OUT}/tmp_sort"

export APPTAINER_HOME="${PWD}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}" "${TMP}"

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${BIND}" "${SIF}" "$@"; }

NC=$(comm -12 \
    <(cut -f1 "${GENOME}.fai" | sort -u) \
    <(awk -F'\t' '$0 !~ /^#/{print $1}' "${TE_GFF}" | sort -u) | wc -l)
echo "[$(date)] contigs communs genome/TE : ${NC}"
[ "${NC}" -eq 0 ] && { echo "ERREUR: TE_GFF sur un autre assemblage." >&2; exit 1; }

run bash -c "export TMPDIR='${TMP}'; set -o pipefail; \
    cut -f1,2 '${GENOME}.fai' > '${TMP}/genome.txt'; \
    sort -T '${TMP}' -k1,1 -k2,2n '${OUT}/mask_mappable.bed' > '${TMP}/mappable.sorted.bed'; \
    awk -F'\t' '\$0 !~ /^#/{print \$1\"\t\"\$4-1\"\t\"\$5}' '${TE_GFF}' \
      | sort -T '${TMP}' -k1,1 -k2,2n > '${TMP}/te.sorted.bed'"

run bash -c "export TMPDIR='${TMP}'; set -o pipefail; \
    bedtools subtract -sorted -g '${TMP}/genome.txt' \
        -a '${TMP}/mappable.sorted.bed' -b '${TMP}/te.sorted.bed' \
    | bedtools merge -i - > '${OUT}/accessible.bed'"

[ -s "${OUT}/accessible.bed" ] || { echo "ERREUR: accessible.bed vide, voir le log." >&2; exit 1; }

# Rapport
TOTAL=$(awk '{s+=$2} END{print s}' "${GENOME}.fai")
KEPT=$(awk '{s+=$3-$2} END{print s}' "${OUT}/accessible.bed")
awk -v k="${KEPT}" -v t="${TOTAL}" 'BEGIN{printf "[final] accessibles: %d / %d pb (%.1f%%)\n", k, t, 100*k/t}'
echo "[$(date)] -> ${OUT}/accessible.bed"