#!/bin/bash
#SBATCH --job-name=hapcaller
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=30:00:00
#SBATCH --output=logs/03_hapcaller_%j.log

set -euo pipefail
SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
INDIR=$(readlink -f results/02_markdup)
OUTDIR="results/03_gvcf"
SAMPLE=${1:?usage: sbatch 03_haplotype_caller.sh <sample_name>}
T=${SLURM_CPUS_PER_TASK:-4}
BIND=$PWD
run()    { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"; }
mkdir -p "${OUTDIR}" logs

# Decompression (once)
REF_REL=$(run yq -r '.reference_genome' "${YAML}")
REF_GZ=$(readlink -f "${REF_REL}")
[[ -s "${REF_GZ}" ]] || { echo "ERROR: reference not found: ${REF_GZ}" >&2; exit 1; }
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

# Index
[[ -s "${REF}.fai" ]] || run samtools faidx "${REF}"
DICT="${REF%.fa}.dict"
[[ -s "${DICT}" ]] || run gatk CreateSequenceDictionary -R "${REF}" -O "${DICT}"

# Check for rerun
BAM="${INDIR}/${SAMPLE}.dedup.bam"
GVCF="${OUTDIR}/${SAMPLE}.g.vcf.gz"
[[ -s "${BAM}" ]] || { echo "ERROR: missing ${BAM}" >&2; exit 1; }
if [[ -s "${GVCF}" ]] && [[ -s "${GVCF}.tbi" ]]; then
    echo "gVCF already exists for ${SAMPLE}, skipping"
    exit 0
fi

# Only keep the ragtag contigs
INTERVAL_FLAGS=""
N_CONTIGS=0
while read -r chrom; do
    INTERVAL_FLAGS="${INTERVAL_FLAGS} -L ${chrom}"
    N_CONTIGS=$((N_CONTIGS+1))
done < <(awk '$1 ~ /^chr.*_RagTag$/ {print $1}' "${REF}.fai")

if [[ "${N_CONTIGS}" -eq 0 ]]; then
    echo "ERROR: no contigs matching ^chr.*_RagTag in ${REF}.fai" >&2
    echo "First 10 contigs:" >&2
    head -10 "${REF}.fai" | cut -f1 >&2
    exit 1
fi
echo "Calling on ${N_CONTIGS} chr*_RagTag scaffolds"


# HaplotypeCaller
# 12G cgroup, -Xmx8G leaves ~4G for JVM off-heap and native pair-HMM.
echo "HaplotypeCaller for ${SAMPLE} at $(date)"
run gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R "${REF}" \
    -I "${BAM}" \
    -O "${GVCF}" \
    -ERC GVCF \
    --native-pair-hmm-threads "${T}" \
    ${INTERVAL_FLAGS}

[[ -s "${GVCF}" && -s "${GVCF}.tbi" ]] || \
    { echo "ERROR: gVCF output incomplete, keeping dedup.bam" >&2; exit 1; }

VARIANTS=$(run_sh "bcftools view -v snps,indels -H '${GVCF}' | wc -l")
echo "Called ${VARIANTS} variant records (SNPs + indels) for ${SAMPLE}"
[[ "${VARIANTS}" -gt 0 ]] || \
    { echo "ERROR: zero variants called, keeping dedup.bam" >&2; exit 1; }

run bcftools stats "${GVCF}" | grep '^SN' | head -20

echo "gVCF OK, removing dedup.bam"
rm -f "${BAM}" "${BAM}.bai"

echo "Done for ${SAMPLE} at $(date)"
echo "Output: ${GVCF}"