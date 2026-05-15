#!/bin/bash

set -euo pipefail

# === Inputs ===
SIF=$(readlink -f images/sif/mitohifi.sif)
ASM=$(readlink -f results/assembly.fa)
MITO_FA=$(readlink -f reference_data/lolium_perenne_mitochondrion.fasta)
MITO_GB=$(readlink -f reference_data/lolium_perenne_mitochondrion.gb)
CHLORO_FA=$(readlink -f reference_data/lolium_perenne_chloroplast.fasta)
CHLORO_GB=$(readlink -f reference_data/lolium_perenne_chloroplast.gb)

# === Output ===
OUTDIR="$(pwd)/results/05_mitohifi"
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)
CACHE="${ROOT}/.cache_mitohifi"

mkdir -p "${CACHE}/home" "${CACHE}/matplotlib" "${OUTDIR}" logs

run() {
    singularity exec \
        --bind "${ROOT}" \
        --env HOME="${CACHE}/home" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

# Run MitoHiFi once for a given organelle.
# Args: label (mito|chloro), reference fasta, reference genbank, NCBI translation code
run_mitohifi() {
    local label=$1 ref_fa=$2 ref_gb=$3 code=$4
    local workdir="${OUTDIR}/${label}_work"

    mkdir -p "${workdir}"
    cd "${workdir}"

    run mitohifi.py \
        -c "${ASM}" \
        -f "${ref_fa}" \
        -g "${ref_gb}" \
        -a plant \
        -o "${code}" \
        -t "${T}"

    # Collect flagged contig IDs
    tail -n+2 contigs_stats.tsv | cut -f1 | sort -u > "${OUTDIR}/${label}.ids"

    # Collect annotated, circularized genome
    cp final_mitogenome.fasta "${OUTDIR}/lmultiflorum_${label}.fa"
    cp final_mitogenome.gb    "${OUTDIR}/lmultiflorum_${label}.gb"

    cd "${ROOT}"
}

echo "=== MitoHiFi started: $(date) ==="

run_mitohifi chloro "${CHLORO_FA}" "${CHLORO_GB}" 11
run_mitohifi mito   "${MITO_FA}"   "${MITO_GB}"   1

# Combined list of organellar contigs
cat "${OUTDIR}/mito.ids" "${OUTDIR}/chloro.ids" | sort -u > "${OUTDIR}/organellar_ids.txt"

echo "=== MitoHiFi done: $(date) ==="
echo "Flagged contigs: $(wc -l < "${OUTDIR}/organellar_ids.txt")"
echo "Mito genome:    ${OUTDIR}/lmultiflorum_mito.fa(.gb)"
echo "Chloro genome:  ${OUTDIR}/lmultiflorum_chloro.fa(.gb)"