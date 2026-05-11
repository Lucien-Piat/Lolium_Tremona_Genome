#!/bin/bash
#SBATCH --job-name=hapcaller
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=36:00:00
#SBATCH --output=logs/03_hapcaller_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
INDIR=$(readlink -f results/02_markdup)
OUTDIR="results/03_gvcf"
SAMPLE=${1:?usage: sbatch 03_haplotypecaller.sh <sample_name>}
T=${SLURM_CPUS_PER_TASK:-4}

BIND=$PWD
run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }

mkdir -p "${OUTDIR}" logs

REF_REL=$(run yq -r '.reference_genome' "${YAML}")
REF_GZ=$(readlink -f "${REF_REL}")
[[ -s "${REF_GZ}" ]] || { echo "ERROR: reference not found: ${REF_GZ}" >&2; exit 1; }

# Reference must be decompressed for GATK (.fa + .fai + .dict)
REF="${REF_GZ%.gz}"
if [[ ! -s "${REF}" ]]; then
    echo "Decompressing reference at $(date)"
    (
        flock -x 200
        if [[ ! -s "${REF}" ]]; then
            run pigz -dcp "${T}" "${REF_GZ}" > "${REF}"
        fi
    ) 200>"${REF}.lock"
    rm -f "${REF}.lock"
fi

[[ -s "${REF}.fai" ]] || run samtools faidx "${REF}"

DICT="${REF%.fa}.dict"
[[ -s "${DICT}" ]] || run gatk CreateSequenceDictionary -R "${REF}" -O "${DICT}"

BAM="${INDIR}/${SAMPLE}.dedup.bam"
GVCF="${OUTDIR}/${SAMPLE}.g.vcf.gz"

[[ -s "${BAM}" ]] || { echo "ERROR: missing ${BAM}" >&2; exit 1; }

if [[ -s "${GVCF}" ]] && [[ -s "${GVCF}.tbi" ]]; then
    echo "gVCF already exists for ${SAMPLE}, skipping"
    exit 0
fi

echo "Running HaplotypeCaller for ${SAMPLE} at $(date)"

# Restrict to chromosome scaffolds (skip unplaced)
INTERVAL_FLAGS=""
for chrom in $(awk '$1 ~ /^chr/ {print $1}' "${REF}.fai"); do
    INTERVAL_FLAGS="${INTERVAL_FLAGS} -L ${chrom}"
done

if [[ -z "${INTERVAL_FLAGS}" ]]; then
    echo "WARN: no contigs matching ^chr in ${REF}.fai, calling on whole genome"
fi

run gatk --java-options "-Xmx40G" HaplotypeCaller \
    -R "${REF}" \
    -I "${BAM}" \
    -O "${GVCF}" \
    -ERC GVCF \
    --native-pair-hmm-threads "${T}" \
    ${INTERVAL_FLAGS}

if [[ ! -s "${GVCF}" ]] || [[ ! -s "${GVCF}.tbi" ]]; then
    echo "ERROR: gVCF output incomplete, keeping dedup.bam" >&2
    exit 1
fi

VARIANTS=$(run bcftools view -v snps,indels --no-version "${GVCF}" | run bcftools view -H | wc -l)
echo "Called ${VARIANTS} variant records for ${SAMPLE}"

if [[ "${VARIANTS}" -le 0 ]]; then
    echo "ERROR: zero variants called, keeping dedup.bam" >&2
    exit 1
fi

# Quick stats
run bcftools stats "${GVCF}" | grep '^SN' | head -20

echo "gVCF OK, removing dedup.bam"
rm -f "${BAM}" "${BAM}.bai"

echo "Done for ${SAMPLE} at $(date)"
echo "Output: ${GVCF}"