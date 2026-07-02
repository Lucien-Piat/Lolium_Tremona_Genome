#!/bin/bash
#SBATCH --job-name=hapcaller
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000M        
#SBATCH --time=05:00:00
#SBATCH --output=logs/03_hapcaller_%j.log

set -euo pipefail
SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
INDIR=$(readlink -f results/02_markdup)
OUTDIR="results/03_gvcf_T"
SAMPLE=${1:?usage: sbatch 03_haplotype_caller.sh <sample_name>}
CPUS=${SLURM_CPUS_PER_TASK:-16}
SCATTER_N=48
BIND=$PWD
run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
mkdir -p "${OUTDIR}" logs

REF_REL=$(run yq -r '.reference_genome' "${YAML}")
REF_GZ=$(readlink -f "${REF_REL}")
REF="${REF_GZ%.gz}"
DICT="${REF%.fa}.dict"
SCATTER="${OUTDIR}/scatter"

(
    flock -x 200
    if [[ ! -s "${REF}" ]]; then
        run pigz -dcp "${CPUS}" "${REF_GZ}" > "${REF}"
    fi
    [[ -s "${REF}.fai" ]] || run samtools faidx "${REF}"
    [[ -s "${DICT}" ]]    || run gatk CreateSequenceDictionary -R "${REF}" -O "${DICT}"
    if [[ ! -s "${SCATTER}/0000-scattered.interval_list" ]]; then
        rm -rf "${SCATTER}"
        run gatk SplitIntervals -R "${REF}" --scatter-count "${SCATTER_N}" -O "${SCATTER}"
    fi
) 200>"${OUTDIR}/.refprep.lock"

BAM="${INDIR}/${SAMPLE}.dedup.bam"
GVCF="${OUTDIR}/${SAMPLE}.g.vcf.gz"
TMP="${OUTDIR}/tmp_${SAMPLE}_${SLURM_JOB_ID:-$$}"
rm -rf "${TMP}"; mkdir -p "${TMP}"
trap 'rm -rf "${TMP}"' EXIT

export SIF BIND REF BAM TMP
printf '%s\0' "${SCATTER}"/*-scattered.interval_list \
| xargs -0 -P "${CPUS}" -I{} bash -c '
    il="$1"
    id=$(basename "${il}" -scattered.interval_list)
    singularity exec --bind "${BIND}" "${SIF}" \
        gatk --java-options "-Xmx2G" HaplotypeCaller \
            -R "${REF}" \
            -I "${BAM}" \
            -O "${TMP}/${id}.g.vcf.gz" \
            -L "${il}" \
            -ERC GVCF \
            --native-pair-hmm-threads 1
' _ {}
# xargs exits non-zero if any interval failed -> set -e stops here

# --- concat shards in genomic order --------------------------------------
ls "${TMP}"/*.g.vcf.gz | sort > "${TMP}/shards.list"
[[ "$(wc -l < "${TMP}/shards.list")" -eq "${SCATTER_N}" ]] \
    || { echo "ERROR: $(wc -l < "${TMP}/shards.list")/${SCATTER_N} shards seulement" >&2; exit 1; }
run bcftools concat -f "${TMP}/shards.list" -Oz -o "${GVCF}"
run bcftools index -t "${GVCF}"

[[ -s "${GVCF}" && -s "${GVCF}.tbi" ]] || \
    { echo "ERROR: gVCF output incomplete, keeping dedup.bam" >&2; exit 1; }

# --- validation (unchanged) ----------------------------------------------
VARIANTS=$(run bcftools view -v snps,indels -H "${GVCF}" | wc -l)
echo "Called ${VARIANTS} variant records (SNPs + indels) for ${SAMPLE}"
[[ "${VARIANTS}" -gt 0 ]] || \
    { echo "ERROR: zero variants called, keeping dedup.bam" >&2; exit 1; }

run bcftools stats "${GVCF}" | grep '^SN' | head -20

rm -f "${BAM}" "${BAM}.bai"
echo "done"