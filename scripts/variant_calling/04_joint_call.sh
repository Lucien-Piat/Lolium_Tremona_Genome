#!/bin/bash
#SBATCH --job-name=joint_call
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/04_joint_call_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
INDIR=$(readlink -f results/03_gvcf)
OUTDIR="results/04_joint_calling"
BIND=$PWD
run()    { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"; }
mkdir -p "${OUTDIR}" logs
REF=$(readlink -f "$(run yq -r '.reference_genome' "${YAML}")")
REF="${REF%.gz}"

[[ -s "${REF}" && -s "${REF}.fai" && -s "${REF%.fa}.dict" ]] \
    || { echo "ERROR: reference, .fai, or .dict missing" >&2; exit 1; }

V_ARGS=""; N=0

for g in "${INDIR}"/*.g.vcf.gz; do
    [[ -s "$g" ]] && { V_ARGS+=" -V $g"; N=$((N+1)); }
done

[[ "${N}" -ge 2 ]] || { echo "ERROR: need >=2 gVCFs, found ${N}" >&2; exit 1; }

echo "Joint calling on ${N} samples"

cd "${OUTDIR}"

[[ -s cohort_combined.g.vcf.gz ]] || \
    run gatk --java-options "-Xmx12G" CombineGVCFs \
        -R "${REF}" ${V_ARGS} -O cohort_combined.g.vcf.gz

[[ -s cohort_raw.vcf.gz ]] || \
    run gatk --java-options "-Xmx12G" GenotypeGVCFs \
        -R "${REF}" -V cohort_combined.g.vcf.gz -O cohort_raw.vcf.gz

# Stritt et al. (2022): SNPs only, QUAL > 20, QD > 8
run_sh "bcftools view -v snps cohort_raw.vcf.gz \
        | bcftools filter -i 'QUAL>20 && INFO/QD>8' -Oz -o cohort_snps_filtered.vcf.gz"
run bcftools index -t cohort_snps_filtered.vcf.gz
run bcftools stats cohort_snps_filtered.vcf.gz > cohort_snps_filtered.stats

RAW=$(run_sh  "bcftools view -H cohort_raw.vcf.gz | wc -l")
SNPS=$(run_sh "bcftools view -H cohort_snps_filtered.vcf.gz | wc -l")
echo "Raw: ${RAW}  SNPs retained: ${SNPS}"