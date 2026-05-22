#!/bin/bash
# Self-synteny (MCScanX) + NUMT/NUPT (minimap2) data prep for Circos.
# Each step runs only if its output file is missing. To redo a step,
# delete its output. To redo everything, rm -rf the OUTDIR.

set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
GENOME="${ROOT}/results/data_circo/lmultiflorum.tremona.placed.fa"
GFF="${ROOT}/results/data_circo/tremona.gene_annotation.placed.gff"
MITO="${ROOT}/results/data_circo/lmul_tremona.mito.fasta"
PLTD="${ROOT}/results/data_circo/lmul_tremona.pltd.fasta"

LIB="${ROOT}/scripts/qc/circo_plot/lib"
DATA_DIR="${ROOT}/results/data_circo"
OUTDIR="${DATA_DIR}/mcscanx"
mkdir -p "${OUTDIR}"

T=4
PREFIX="lm"
NAME="tremona"

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

cd "${OUTDIR}"

# Stage 1: proteins + MCScanX gff
if [ ! -s proteins.raw.faa ]; then
    echo "[$(date)] gffread"
    run gffread -y proteins.raw.faa -g "${GENOME}" "${GFF}"
fi

if [ ! -s proteins.faa ]; then
    echo "[$(date)] clean proteins"
    run python3 "${LIB}/clean_proteins.py" proteins.raw.faa proteins.faa
fi

if [ ! -s "${NAME}.gff" ]; then
    echo "[$(date)] MCScanX gff"
    run awk -v p="${PREFIX}" -F'\t' '
        $3 == "mRNA" {
            match($9, /ID=([^;]+)/, a)
            chrom = $1; sub(/^chr/, p, chrom)
            print chrom"\t"a[1]"\t"$4"\t"$5
        }' "${GFF}" > "${NAME}.gff"
fi

# Stage 2: DIAMOND
if [ ! -s proteins.db.dmnd ]; then
    echo "[$(date)] diamond makedb"
    run diamond makedb --in proteins.faa -d proteins.db --threads "${T}"
fi

if [ ! -s "${NAME}.blast" ]; then
    echo "[$(date)] diamond blastp"
    run diamond blastp \
        -q proteins.faa -d proteins.db -o "${NAME}.blast" \
        -e 1e-10 --outfmt 6 --max-target-seqs 5 --more-sensitive \
        --block-size 1.5 --index-chunks 2 --threads "${T}"
fi

# Stage 3: MCScanX + links
if [ ! -s "${NAME}.collinearity" ]; then
    echo "[$(date)] MCScanX"
    run MCScanX -s 5 -e 1e-10 -m 25 "${OUTDIR}/${NAME}"
fi

if [ ! -s self_synteny_links.tsv ]; then
    echo "[$(date)] collinearity to links"
    run python3 "${LIB}/collinearity_to_links.py" \
        "${NAME}" "${PREFIX}" self_synteny_links.tsv
fi

# Stage 4: NUMT / NUPT
cd "${DATA_DIR}"

map_organelle() {
    local org="$1" tag="$2" low="${2,,}"
    local paf="${low}_to_nuclear.paf"
    local links="${low}_links.tsv"

    if [ ! -s "${paf}" ]; then
        echo "[$(date)] minimap2 ${tag}"
        run minimap2 -t "${T}" -cx asm20 -N 50 --secondary=yes -K 10M \
            "${GENOME}" "${org}" > "${paf}"
    fi

    if [ ! -s "${links}" ]; then
        echo "[$(date)] filter ${tag}"
        run awk -v ml=500 -v mi=0.80 -v tag="${tag}" '
            $11 >= ml && ($10/$11) >= mi {
                print $1, $3, $4, $6, $8, $9, $5, tag, ($10/$11)
            }' OFS='\t' "${paf}" > "${links}"
    fi
}

map_organelle "${MITO}" "NUMT"
map_organelle "${PLTD}" "NUPT"

if [ ! -s organelle_insertion_links.tsv ]; then
    cat numt_links.tsv nupt_links.tsv > organelle_insertion_links.tsv
fi

# Stats (always rewritten)
{
    echo "===== self synteny ====="
    echo "Blocks: $(tail -n +2 ${OUTDIR}/self_synteny_links.tsv | wc -l)"
    tail -n +2 "${OUTDIR}/self_synteny_links.tsv" \
        | awk '{print $1"\t"$4}' | sort | uniq -c | sort -rn | head -10
    echo
    echo "===== NUMTs ====="; wc -l < numt_links.tsv
    echo "===== NUPTs ====="; wc -l < nupt_links.tsv
} | tee synteny_stats.txt

echo "[$(date)] Done"