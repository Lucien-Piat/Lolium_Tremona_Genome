#!/bin/bash
#SBATCH --job-name=mitohifi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --output=logs/05_mitohifi_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/polish.sif)
ASM_GZ=$(readlink -f results/03_purge/lmultiflorum.purged.fa.gz)
OUTDIR="results/04_organellar"
TRACKING="results/assembly_tracking.tsv"
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)

BIND="/cluster/scratch"

MITO_FA=$(readlink -f reference_data/lolium_perenne_mitochondrion.fasta)
MITO_GB=$(readlink -f reference_data/lolium_perenne_mitochondrion.gb)
CHLORO_FA=$(readlink -f reference_data/lolium_perenne_chloroplast.fasta)
CHLORO_GB=$(readlink -f reference_data/lolium_perenne_chloroplast.gb)

CACHE="$(pwd)/.cache_mitohifi"
mkdir -p "${CACHE}/home" "${CACHE}/matplotlib" "${CACHE}/fontconfig" "${OUTDIR}" logs

run() {
    singularity exec \
        --bind "${BIND}" \
        --env HOME="${CACHE}/home" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env FONTCONFIG_PATH="${CACHE}/fontconfig" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

track() {
    local stage=$1 fa=$2
    local stats nseq size
    stats=$(run seqkit stats -T "${fa}" | tail -1)
    nseq=$(echo "${stats}" | cut -f4)
    size=$(echo "${stats}" | cut -f5)
    printf '%s\t%s\t%s\t%s\n' "${stage}" "$(readlink -f "${fa}")" "${nseq}" "${size}" >> "${TRACKING}"
    echo "${stage}: ${nseq} contigs, ${size} bp"
}

run_mitohifi() {
    local label=$1 ref_fa=$2 ref_gb=$3 code=$4
    local workdir="${ROOT}/${OUTDIR}/${label}"
    local ids_file="${ROOT}/${OUTDIR}/${label}.ids"
    local asm_abs="${ROOT}/${OUTDIR}/assembly.fa"

    mkdir -p "${workdir}"
    > "${ids_file}"
    cd "${workdir}"

    if run mitohifi.py \
            -c "${asm_abs}" \
            -f "${ref_fa}" \
            -g "${ref_gb}" \
            -o "${code}" -t "${T}"; then
        echo "INFO: mitohifi ${label} run completed"
    else
        echo "WARN: mitohifi ${label} exited non-zero, checking outputs anyway"
    fi

    if [[ -f contigs_stats.tsv ]] && [[ $(tail -n+2 contigs_stats.tsv | wc -l) -gt 0 ]]; then
        tail -n+2 contigs_stats.tsv | cut -f1 | sort -u > "${ids_file}"
        [[ -f final_mitogenome.fasta ]] && cp final_mitogenome.fasta "${ROOT}/${OUTDIR}/lmultiflorum_${label}.fa"
        [[ -f final_mitogenome.gb ]]    && cp final_mitogenome.gb    "${ROOT}/${OUTDIR}/lmultiflorum_${label}.gb"
        echo "INFO: ${label} contigs identified: $(wc -l < "${ids_file}")"
    else
        echo "WARN: no ${label} contigs identified"
    fi

    cd "${ROOT}"
    rm -rf "${workdir}"
}

[[ -f "${TRACKING}" ]] || printf 'stage\tfile\tcontigs\tsize\n' > "${TRACKING}"

run pigz -dcp "${T}" "${ASM_GZ}" > "${OUTDIR}/assembly.fa"

run_mitohifi mito   "${MITO_FA}"   "${MITO_GB}"   1
run_mitohifi chloro "${CHLORO_FA}" "${CHLORO_GB}" 11

cd "${OUTDIR}"
cat mito.ids chloro.ids 2>/dev/null | sort -u > organellar_ids.txt

if [[ -s organellar_ids.txt ]]; then
    run seqkit grep -v -f organellar_ids.txt assembly.fa > lmultiflorum.nuclear.fa
    echo "INFO: $(wc -l < organellar_ids.txt) organellar contigs removed"
else
    echo "WARN: no organellar contigs found, keeping assembly as-is"
    mv assembly.fa lmultiflorum.nuclear.fa
fi

rm -f assembly.fa mito.ids chloro.ids organellar_ids.txt

for f in lmultiflorum_mito.fa lmultiflorum_chloro.fa lmultiflorum.nuclear.fa; do
    [[ -f "${f}" ]] && run pigz -p "${T}" "${f}"
done

cd "${ROOT}"
rm -rf "${CACHE}"
track "assembly-nuclear" "${OUTDIR}/lmultiflorum.nuclear.fa.gz"