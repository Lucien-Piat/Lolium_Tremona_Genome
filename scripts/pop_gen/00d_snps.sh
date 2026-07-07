#!/bin/bash
#SBATCH --job-name=filt_snps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --time=04:00:00
#SBATCH --output=logs/filter_%x_%j.log
set -euo pipefail
# SNP : masque + bialleliques + MAF, puis pruning LD (plink2), cible ~50k.

ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
SNPS="${ROOT}/results/04_joint_calling/cohort_snps_filtered.vcf.gz"
MASK="${ROOT}/results/mask/accessible.bed"
OUT="${ROOT}/results/filtered_vcf"

MAF=0.05         # a confirmer avec la cheffe, pas 0.5
LD_R2=0.1        # seuil de pruning LD (retire les SNP correles au-dessus)
LD_WIN=50        # fenetre de pruning en kb
LD_STEP=1        # pas obligatoirement 1 quand la fenetre est en kb
T="${SLURM_CPUS_PER_TASK:-4}"
MEM_MB=$(( ${SLURM_MEM_PER_CPU:-1024} * T - 1024 ))   # headroom de 1 Go pour plink2

export APPTAINER_HOME="${ROOT}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}" "${OUT}" logs

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

#########################################################
# Etape 1 : masque + bialleliques + MAF
#########################################################
if [ ! -s "${OUT}/snps.masked.biallelic.vcf.gz" ]; then
    echo "[$(date)] SNP : masque + bialleliques + MAF ${MAF}"
    run bcftools view -R "${MASK}" -m2 -M2 -v snps \
        -i "MAF>=${MAF}" --threads "${T}" \
        "${SNPS}" -Oz -o "${OUT}/snps.masked.biallelic.vcf.gz"
    run tabix -p vcf "${OUT}/snps.masked.biallelic.vcf.gz"
fi
N_BEFORE=$(run bash -c "bcftools index -n '${OUT}/snps.masked.biallelic.vcf.gz'")
echo "[$(date)] SNP apres MAF : ${N_BEFORE}"

#########################################################
# Etape 2 : pruning LD (plink2)
#########################################################
# Fenetre en kb -> le pas DOIT valoir 1 (contrainte plink2).
if [ ! -s "${OUT}/snps.ld.vcf.gz" ]; then
    echo "[$(date)] SNP : pruning LD (r2>${LD_R2}, fenetre ${LD_WIN}kb pas ${LD_STEP})"
    run plink2 --vcf "${OUT}/snps.masked.biallelic.vcf.gz" \
        --allow-extra-chr --set-all-var-ids '@:#' \
        --threads "${T}" --memory "${MEM_MB}" \
        --indep-pairwise "${LD_WIN}kb" "${LD_STEP}" "${LD_R2}" \
        --out "${OUT}/ldprune"

    run plink2 --vcf "${OUT}/snps.masked.biallelic.vcf.gz" \
        --allow-extra-chr --set-all-var-ids '@:#' \
        --threads "${T}" --memory "${MEM_MB}" \
        --extract "${OUT}/ldprune.prune.in" \
        --export vcf bgz --out "${OUT}/snps.ld"
    run tabix -p vcf "${OUT}/snps.ld.vcf.gz"
fi
N_AFTER=$(run bash -c "bcftools index -n '${OUT}/snps.ld.vcf.gz'")
echo "[$(date)] SNP apres pruning LD : ${N_AFTER}  (cible ~50000)"

echo "[$(date)] Termine. -> ${OUT}/snps.ld.vcf.gz"