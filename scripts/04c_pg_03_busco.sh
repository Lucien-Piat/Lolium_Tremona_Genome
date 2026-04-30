#!/bin/bash
#SBATCH --job-name=pg_busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/04c_pg_03_busco_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/purgegrass.sif)
OUTDIR="results/04c_purgegrass"
DB_ARCHIVE=$(readlink -f reference_data/busco_downloads.tar.gz)
LINEAGE="poales_odb12"
T=${SLURM_CPUS_PER_TASK:-4}

BIND="/cluster/scratch"

CACHE="$(pwd)/.cache_pg_busco"
mkdir -p "${CACHE}/home" "${CACHE}/matplotlib"

run() {
    singularity exec --bind "${BIND}" \
        --env HOME="${CACHE}/home" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

cd "${OUTDIR}"

TMPDB=".tmpdb_busco"
if [[ ! -d "${TMPDB}/lineages/${LINEAGE}" ]]; then
    rm -rf "${TMPDB}"
    mkdir -p "${TMPDB}"
    echo "Extracting BUSCO DB..."
    run tar xzf "${DB_ARCHIVE}" -C "${TMPDB}"
fi

# Verify the extracted DB has the lineage we asked for
if [[ ! -d "${TMPDB}/lineages/${LINEAGE}" ]]; then
    echo "ERROR: lineage ${LINEAGE} not found in ${DB_ARCHIVE}" >&2
    echo "Available lineages:" >&2
    ls "${TMPDB}/lineages/" >&2
    exit 1
fi

# Adjust the lineage path if the archive structure has busco_downloads/ as a top dir
LINEAGE_PATH="${TMPDB}/lineages/${LINEAGE}"
if [[ ! -d "${LINEAGE_PATH}" ]] && [[ -d "${TMPDB}/busco_downloads/lineages/${LINEAGE}" ]]; then
    LINEAGE_PATH="${TMPDB}/busco_downloads/lineages/${LINEAGE}"
fi

echo "Running BUSCO at $(date) with ${T} threads"
run busco \
    -i assembly.fa \
    -o busco_out \
    -l "${LINEAGE_PATH}" \
    -m genome \
    -c "${T}" \
    --offline \
    -f

echo "BUSCO done at $(date)"
echo "Summary:"
grep -E '^\s+C:' busco_out/short_summary*.txt | head -1
echo ""
echo "full_table.tsv path needed for PurgeGrass:"
find busco_out -name 'full_table.tsv' | head -1

rm -rf "${TMPDB}" "${CACHE}"