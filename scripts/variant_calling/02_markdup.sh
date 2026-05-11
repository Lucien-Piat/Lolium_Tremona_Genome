#!/bin/bash
#SBATCH --job-name=markdup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2500M
#SBATCH --time=6:00:00
#SBATCH --output=logs/02_markdup_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/varcall.sif)
INDIR=$(readlink -f results/01_aligned)
OUTDIR="results/02_markdup"
SAMPLE=${1:?usage: sbatch 02_markdup.sh <sample_name>}
T=${SLURM_CPUS_PER_TASK}
BIND="/cluster/scratch"

run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }

mkdir -p "${OUTDIR}" logs
TMP="${OUTDIR}/tmp_${SAMPLE}"
mkdir -p "${TMP}"

SORTED="${INDIR}/${SAMPLE}.sorted.bam"

cd "${OUTDIR}"

# Mark duplicates
run sambamba markdup \
    -t "${T}" \
    --tmpdir="${TMP}" \
    --hash-table-size=4194304 \
    --overflow-list-size=1000000 \
    --sort-buffer-size=16384 \
    --io-buffer-size=512 \
    "${SORTED}" \
    "${SAMPLE}.dedup.bam"

run samtools flagstat -@ "${T}" "${SAMPLE}.dedup.bam" > "${SAMPLE}.dedup.flagstat.txt"

# Verify dedup BAM is sane before cleaning up
if [[ ! -s "${SAMPLE}.dedup.bam" ]] || [[ ! -s "${SAMPLE}.dedup.bam.bai" ]]; then
    echo "ERROR: dedup output incomplete, keeping inputs" >&2
    exit 1
fi

ALIGNED=$(grep "in total" "${SAMPLE}.dedup.flagstat.txt" | head -1 | awk '{print $1}')
if [[ -z "${ALIGNED}" ]] || [[ "${ALIGNED}" -le 0 ]]; then
    echo "ERROR: dedup BAM has no aligned reads, keeping inputs" >&2
    exit 1
fi

echo "dedup OK (${ALIGNED} reads), removing sorted.bam"
rm -f "${SORTED}" "${SORTED}.bai"

rm -rf "${TMP}"
echo "done"