#!/bin/bash
#SBATCH --job-name=mm2_align
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1250M
#SBATCH --time=08:00:00
#SBATCH --output=logs/01_align_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
OUTDIR="results/01_aligned"
SAMPLE=${1:?usage: sbatch 01_align.sh <sample_name>}
T=${SLURM_CPUS_PER_TASK:-4}
BIND=$PWD

run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }

mkdir -p "${OUTDIR}" logs

REF_REL=$(run yq -r '.reference_genome' "${YAML}")
REF=$(readlink -f "${REF_REL}")
MMI="${REF%.fa.gz}.sr.mmi"

[[ -s "${REF}" ]] || { echo "ERROR: reference not found: ${REF}" >&2; exit 1; }

if [[ ! -s "${MMI}" ]]; then
    echo "Index not found, building ${MMI} at $(date)"
    (
        flock -x 200
        if [[ ! -s "${MMI}" ]]; then
            run minimap2 -t "${T}" -x sr -d "${MMI}" "${REF}"
            echo "Index built at $(date)"
        else
            echo "Index built by another job while we were waiting"
        fi
    ) 200>"${MMI}.lock"
    rm -f "${MMI}.lock"
else
    echo "Reusing existing index ${MMI}"
fi

R1=$(readlink -f "$(run yq -r ".samples[] | select(.name == \"${SAMPLE}\") | .r1" "${YAML}")")
R2=$(readlink -f "$(run yq -r ".samples[] | select(.name == \"${SAMPLE}\") | .r2" "${YAML}")")

cd "${OUTDIR}"

# Align, sort, index
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
run bash -c "minimap2 -ax sr -t $((T-2)) -R '${RG}' '${MMI}' '${R1}' '${R2}' \
    | samtools sort -@ 2 -m 2G -O bam -o '${SAMPLE}.sorted.bam' -"
run samtools index -@ "${T}" "${SAMPLE}.sorted.bam"
run samtools flagstat -@ "${T}" "${SAMPLE}.sorted.bam" > "${SAMPLE}.flagstat.txt"

echo "done"