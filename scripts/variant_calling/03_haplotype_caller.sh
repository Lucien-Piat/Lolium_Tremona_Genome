#!/bin/bash
#SBATCH --job-name=hapcaller
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1500M
#SBATCH --time=50:00:00
#SBATCH --output=logs/03_hapcaller_%j.log

# Carefull the mem / cpu is finetuned for Lmul Tremona (to the gigabate)

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
REF="${REF_GZ%.gz}"
if [[ ! -s "${REF}" ]]; then
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

# Gatk 
run gatk --java-options "-Xmx2G" HaplotypeCaller \
    -R "${REF}" \
    -I "${BAM}" \
    -O "${GVCF}" \
    -ERC GVCF \
    --native-pair-hmm-threads "${T}"

[[ -s "${GVCF}" && -s "${GVCF}.tbi" ]] || \
    { echo "ERROR: gVCF output incomplete, keeping dedup.bam" >&2; exit 1; }

VARIANTS=$(run_sh "bcftools view -v snps,indels -H '${GVCF}' | wc -l")
echo "Called ${VARIANTS} variant records (SNPs + indels) for ${SAMPLE}"
[[ "${VARIANTS}" -gt 0 ]] || \
    { echo "ERROR: zero variants called, keeping dedup.bam" >&2; exit 1; }

run bcftools stats "${GVCF}" | grep '^SN' | head -20

rm -f "${BAM}" "${BAM}.bai"

echo "done"