#!/bin/bash
#SBATCH --job-name=joint_call
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=120:00:00
#SBATCH --output=logs/04_joint_%j.log

set -euo pipefail
SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
GVCFDIR=$(readlink -f results/03_gvcf)
OUTDIR="results/04_joint_calling"
MAP="${OUTDIR}/sample_map.tsv"
T=${SLURM_CPUS_PER_TASK:-1}
BIND=$PWD
run()    { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"; }
mkdir -p "${OUTDIR}" logs
WORKTMP="${OUTDIR}/tmp_${SLURM_JOB_ID:-$$}"
mkdir -p "${WORKTMP}"
trap 'rm -rf "${WORKTMP}"' EXIT
DB="${WORKTMP}/genomicsdb"

REF=$(readlink -f "$(run yq -r '.reference_genome' "${YAML}")")
[[ -s "${REF}" && -s "${REF}.fai" && -s "${REF%.fa}.dict" ]] \
    || { echo "ERROR: reference, .fai ou .dict manquant" >&2; exit 1; }

if [[ ! -s "${MAP}" ]]; then
    for g in "${GVCFDIR}"/*.g.vcf.gz; do
        s=$(basename "${g}" .g.vcf.gz)
        printf "%s\t%s\n" "${s}" "$(readlink -f "${g}")"
    done > "${MAP}"
fi
cut -f1 "${REF}.fai" > "${OUTDIR}/intervals.list"
echo "Cohorte : $(wc -l < "${MAP}") echantillons, $(wc -l < "${OUTDIR}/intervals.list") chromosomes"

# GenomicsDB 
run gatk --java-options "-Xmx12G" GenomicsDBImport \
    --sample-name-map "${MAP}" \
    --genomicsdb-workspace-path "${DB}" \
    -L "${OUTDIR}/intervals.list" \
    --merge-input-intervals \
    --batch-size 50 \
    --reader-threads "${T}" \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir "${WORKTMP}"

# 2. all-sites
ALLSITES="${OUTDIR}/cohort_allsites.vcf.gz"
run gatk --java-options "-Xmx18G" GenotypeGVCFs \
    -R "${REF}" \
    -V "gendb://${DB}" \
    --include-non-variant-sites \
    -O "${ALLSITES}" \
    --tmp-dir "${WORKTMP}"
[[ -s "${ALLSITES}" && -s "${ALLSITES}.tbi" ]] \
    || { echo "ERROR: all-sites incomplet" >&2; exit 1; }

# 3. variants-only filtre : SNPs, QUAL>20, QD>8 (Stritt et al. 2022)
SNPS="${OUTDIR}/cohort_snps_filtered.vcf.gz"
run_sh "bcftools view -v snps '${ALLSITES}' \
        | bcftools filter -i 'QUAL>20 && INFO/QD>8' -Oz -o '${SNPS}'"
run bcftools index -t "${SNPS}"
run bcftools stats "${SNPS}" > "${SNPS%.vcf.gz}.stats"

echo "all-sites       : ${ALLSITES} ($(run bcftools index -n "${ALLSITES}") sites)"
echo "variants filtre : ${SNPS} ($(run bcftools index -n "${SNPS}") SNPs)"