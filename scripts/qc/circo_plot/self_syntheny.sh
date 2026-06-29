#!/bin/bash
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
GENOME="${ROOT}/reference_data/lmultiflorum.tremona.fa"
MITO="${ROOT}/results/data_circo/lmul_tremona.mito.fasta"
PLTD="${ROOT}/results/data_circo/lmul_tremona.pltd.fasta"

DATA_DIR="${ROOT}/results/data_circo"
mkdir -p "${DATA_DIR}"

T=4

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

cd "${DATA_DIR}"

map_mito() {
    local org="$1" tag="$2" minlen="$3" minid="$4"
    local low="${tag,,}"
    local paf="${low}_to_nuclear.paf"
    local links="${low}_links.tsv"

    if [ ! -s "${paf}" ]; then
        echo "[$(date)] minimap2 ${tag}"
        run minimap2 -t "${T}" -cx asm20 -N 50 --secondary=yes -K 10M \
            "${GENOME}" "${org}" > "${paf}"
    fi
    if [ ! -s "${links}" ]; then
        echo "[$(date)] filter ${tag} (minlen=${minlen}, minid=${minid})"
        run awk -v ml="${minlen}" -v mi="${minid}" -v tag="${tag}" '
            $11 >= ml && ($10/$11) >= mi {
                print $1, $3, $4, $6, $8, $9, $5, tag, ($10/$11)
            }' OFS='\t' "${paf}" > "${links}"
    fi
}

map_pldt() {
    local org="$1" tag="$2" minlen="$3" minid="$4"
    local low="${tag,,}"
    local paf="${low}_to_nuclear.paf"
    local links="${low}_links.tsv"

    if [ ! -s "${paf}" ]; then
        echo "[$(date)] minimap2 ${tag}"
        run minimap2 -t "${T}" -k 15 -w 10 -A 1 -B 2 -O 2,32 -E 1,0 \
            -N 50 --secondary=yes -K 10M \
            "${GENOME}" "${org}" > "${paf}"
    fi
    if [ ! -s "${links}" ]; then
        echo "[$(date)] filter ${tag} (minlen=${minlen}, minid=${minid})"
        run awk -v ml="${minlen}" -v mi="${minid}" -v tag="${tag}" '
            $11 >= ml && ($10/$11) >= mi {
                print $1, $3, $4, $6, $8, $9, $5, tag, ($10/$11)
            }' OFS='\t' "${paf}" > "${links}"
    fi
}

map_mito "${MITO}" "NUMT" 500 0.80
map_pldt "${PLTD}" "NUPT" 200 0.70

if [ ! -s organelle_insertion_links.tsv ]; then
    cat numt_links.tsv nupt_links.tsv > organelle_insertion_links.tsv
fi

{
    echo "===== NUMTs ====="; wc -l < numt_links.tsv
    echo "===== NUPTs ====="; wc -l < nupt_links.tsv
} | tee organelle_stats.txt

echo "[$(date)] Done"