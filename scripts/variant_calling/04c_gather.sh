#!/bin/bash
#SBATCH --job-name=jc_gather
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/04_gather_%j.log

# Step 3: concat shards -> all-sites, filter -> SNPs, verify, free shards.
# Run AFTER the shard array: N=200 sbatch scripts/04c_gather.sh
set -euo pipefail
cd "$(readlink -f .)"
N=${N:-200}
SIF=$(readlink -f images/sif/varcall.sif)
OUTDIR="results/04_joint_calling"
run()    { singularity exec --bind "$PWD" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "$PWD" "${SIF}" bash -c "$@"; }

got=$(ls "${OUTDIR}"/shards/*/allsites.vcf.gz 2>/dev/null | wc -l)
[[ "${got}" -eq "${N}" ]] || { echo "ERROR: ${got}/${N} shards seulement" >&2; exit 1; }

ls "${OUTDIR}"/shards/*/allsites.vcf.gz | sort > "${OUTDIR}/vcf.list"
rm -f "${OUTDIR}/cohort_allsites.vcf.gz" "${OUTDIR}/cohort_allsites.vcf.gz.tbi"

ALLSITES="${OUTDIR}/cohort_allsites.vcf.gz"
if ! run_sh "bcftools concat --naive-force -f '${OUTDIR}/vcf.list' -o '${ALLSITES}'"; then
    echo "naive-force indisponible, recompression (--threads 2)" >&2
    run_sh "bcftools concat --threads 2 -f '${OUTDIR}/vcf.list' -Oz -o '${ALLSITES}'"
fi
run bcftools index -t "${ALLSITES}"

N_ALL=$(run bcftools index -s "${ALLSITES}" | awk '{s+=$3} END{print s+0}')
echo "all-sites       : ${ALLSITES} (${N_ALL} sites)"

# variants-only : SNPs, QUAL>20, QD>8 (Stritt et al. 2022)
SNPS="${OUTDIR}/cohort_snps_filtered.vcf.gz"
run_sh "bcftools view -v snps '${ALLSITES}' \
        | bcftools filter -i 'QUAL>20 && INFO/QD>8' --threads 2 -Oz -o '${SNPS}'"
run bcftools index -t "${SNPS}"
N_SNP=$(run bcftools index -n "${SNPS}")
run bcftools stats "${SNPS}" > "${SNPS%.vcf.gz}.stats"
echo "variants filtre : ${SNPS} (${N_SNP} SNPs)"

if [[ "${N_ALL}" -gt 0 && "${N_SNP}" -gt 0 && -s "${ALLSITES}.tbi" && -s "${SNPS}.tbi" ]]; then
    du -sh "${OUTDIR}/shards" 2>/dev/null | awk '{print "shards liberes : "$1}'
    rm -rf "${OUTDIR}/shards"
    rm -f "${OUTDIR}/vcf.list"
    echo "shards supprimes (verif OK)"
else
    echo "WARN: verification incomplete, shards CONSERVES" >&2
    exit 1
fi
echo "done"