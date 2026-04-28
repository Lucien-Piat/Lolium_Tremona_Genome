#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --output=logs/09_busco_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/evaluation.sif)
OUTDIR="results/08_busco"
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)

BIND="/cluster/scratch"
LINEAGE="embryophyta_odb12"

DB_ARCHIVE=$(readlink -f reference_data/busco_downloads.tar.gz)

# Writable caches, $HOME is read-only on compute nodes.
CACHE="$(pwd)/.cache_busco"
mkdir -p "${CACHE}/home" "${CACHE}/matplotlib" "${CACHE}/fontconfig"

run() {
    singularity exec \
        --bind "${BIND}" \
        --env HOME="${CACHE}/home" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env FONTCONFIG_PATH="${CACHE}/fontconfig" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

ASSEMBLIES=(
    primary:results/02_assembly/lmultiflorum.bp.p_ctg.fa.gz
    purged:results/03_purge/lmultiflorum.purged.fa.gz
    nuclear:results/04_organellar/lmultiflorum.nuclear.fa.gz
    decontam:results/05_blobtoolkit/lmultiflorum.decontam.fa.gz
    scaffolded:results/06_scaffolding/lmultiflorum.scaffolded.fa.gz
)

mkdir -p "${OUTDIR}" logs

# Extract the lineage DB once, reuse for all samples.
TMPDB="${OUTDIR}/.tmpdb"
echo "Untarring BUSCO lineage DB..."
mkdir -p "${TMPDB}"
run bash -c "pigz -dc -p ${T} ${DB_ARCHIVE} | tar -xf - -C ${TMPDB}"

for entry in "${ASSEMBLIES[@]}"; do
    label="${entry%%:*}"
    gz="${entry#*:}"
    if [[ ! -f "${gz}" ]]; then
        echo "SKIP: ${label} at ${gz} does not exist"
        continue
    fi

    summary="${OUTDIR}/${label}/short_summary.specific.${LINEAGE}.${label}.txt"
    if [[ -f "${summary}" ]]; then
        echo "SKIP: BUSCO already done for ${label}"
        continue
    fi

    echo "=========================================="
    echo "BUSCO on ${label}"
    echo "=========================================="

    tmpfasta="${OUTDIR}/${label}.fa"
    run pigz -dcp "${T}" "${gz}" > "${tmpfasta}"

    run busco \
        -i "${tmpfasta}" \
        -o "${label}" \
        -l "${TMPDB}/lineages/${LINEAGE}" \
        -m genome \
        -c "${T}" \
        --out_path "${OUTDIR}" \
        --offline \
        -f

    # Keep only the short summary, drop heavy intermediates.
    rm -f "${tmpfasta}"
    find "${OUTDIR}/${label}" -mindepth 1 \
        ! -name 'short_summary.specific.*' -delete 2>/dev/null || true
done

# Clean up the extracted DB and caches.
rm -rf "${TMPDB}" "${CACHE}"

# Compact summary.
echo ""
echo "=============================="
echo "BUSCO summary across assemblies"
echo "=============================="
for summary in "${OUTDIR}"/*/short_summary.specific.*.txt; do
    [[ -f "${summary}" ]] || continue
    label=$(basename "$(dirname "${summary}")")
    line=$(grep -E '^\s+C:' "${summary}" | head -1 | xargs)
    printf '%-15s %s\n' "${label}" "${line}"
done