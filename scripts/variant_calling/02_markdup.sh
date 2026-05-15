#!/bin/bash
#SBATCH --job-name=markdup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=2:00:00
#SBATCH --output=logs/02_markdup_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/varcall.sif)
INDIR=$(readlink -f results/01_aligned)
OUTDIR=$(readlink -f results/02_markdup)
SAMPLE=${1:?usage: sbatch 02_markdup.sh <sample_name>}
T=${SLURM_CPUS_PER_TASK}
BIND="/cluster/scratch"
run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}" logs
TMP="${OUTDIR}/tmp_${SAMPLE}"
mkdir -p "${TMP}"
SORTED="${INDIR}/${SAMPLE}.sorted.bam"
DEDUP="${OUTDIR}/${SAMPLE}.dedup.bam"
FLAGSTAT="${OUTDIR}/${SAMPLE}.dedup.flagstat.txt"

[[ -s "${SORTED}" ]] || { echo "ERROR: missing input ${SORTED}" >&2; exit 1; }

# Mark duplicates
run sambamba markdup \
    -t "${T}" \
    --tmpdir="${TMP}" \
    --hash-table-size=4194304 \
    --overflow-list-size=1000000 \
    --sort-buffer-size=16384 \
    --io-buffer-size=512 \
    "${SORTED}" \
    "${DEDUP}"

run samtools flagstat -@ "${T}" "${DEDUP}" > "${FLAGSTAT}"

# Cleanup
if [[ ! -s "${DEDUP}" ]] || [[ ! -s "${DEDUP}.bai" ]]; then
    echo "ERROR: dedup output incomplete, keeping inputs" >&2
    exit 1
fi

ALIGNED=$(grep "in total" "${FLAGSTAT}" | head -1 | awk '{print $1}')
if [[ -z "${ALIGNED}" ]] || [[ "${ALIGNED}" -le 0 ]]; then
    echo "ERROR: dedup BAM has no aligned reads, keeping inputs" >&2
    exit 1
fi

echo "dedup OK (${ALIGNED} reads), removing sorted.bam"
rm -f "${SORTED}" "${SORTED}.bai"
rm -rf "${TMP}"
echo "done"