#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --output=logs/08_quast_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/evaluation.sif)
OUTDIR="results/07_quast"
T=${SLURM_CPUS_PER_TASK:-4}

BIND="/cluster/scratch"

run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }

# All assemblies QUAST will compare, in pipeline order.
ASSEMBLIES=(
    results/02_assembly/lmultiflorum.bp.p_ctg.fa.gz
    results/02_assembly/lmultiflorum.bp.hap1.p_ctg.fa.gz
    results/02_assembly/lmultiflorum.bp.hap2.p_ctg.fa.gz
    results/03_purge/lmultiflorum.purged.fa.gz
    results/04_organellar/lmultiflorum.nuclear.fa.gz
    results/05_blobtoolkit/lmultiflorum.decontam.fa.gz
    results/06_scaffolding/lmultiflorum.scaffolded.fa.gz
)

mkdir -p "${OUTDIR}" logs

INPUTS=()
LABELS_A=()
for gz in "${ASSEMBLIES[@]}"; do
    if [[ -f "${gz}" ]]; then
        INPUTS+=("${gz}")
        LABELS_A+=("$(basename "${gz}" .fa.gz)")
    else
        echo "SKIP: ${gz} does not exist"
    fi
done
LABELS=$(IFS=, ; echo "${LABELS_A[*]}")

[[ ${#INPUTS[@]} -gt 0 ]] || { echo "ERROR: no assemblies found" >&2; exit 1; }

echo "Running QUAST on ${#INPUTS[@]} assemblies: ${LABELS}"

run quast \
    "${INPUTS[@]}" \
    -o "${OUTDIR}" \
    --threads "${T}" \
    --large \
    --no-snps \
    --no-plots \
    --scaffold-gap-max-size 10 \
    --split-scaffolds \
    --labels "${LABELS}"

rm -f "${OUTDIR}/quast.log" "${OUTDIR}"/*.tex
rm -rf "${OUTDIR}/basic_stats"

echo ""
echo "QUAST report at: ${OUTDIR}/report.html"
cat "${OUTDIR}/report.txt"