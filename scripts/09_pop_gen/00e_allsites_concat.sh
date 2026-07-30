#!/bin/bash
#SBATCH --job-name=concat_as
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=2G
#SBATCH --time=03:00:00
#SBATCH --output=logs/filter_%x_%j.log
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
FAI="${ROOT}/reference_data/lmultiflorum.tremona.fa.fai"
OUT="${ROOT}/results/filtered_vcf"
TMP="${OUT}/tmp_allsites"
T="${SLURM_CPUS_PER_TASK:-4}"

export APPTAINER_HOME="${ROOT}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}"

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

MISSING=0
while read -r chr _; do
    [ -s "${TMP}/${chr}.hc.vcf.gz" ] || { echo "MANQUE: ${chr}" >&2; MISSING=1; }
done < "${FAI}"
[ "${MISSING}" -eq 0 ] || { echo "ERREUR: morceaux manquants, concat annulee." >&2; exit 1; }

cut -f1 "${FAI}" | sed "s|^|${TMP}/|; s|\$|.hc.vcf.gz|" > "${TMP}/concat_list.txt"

run bash -c "set -o pipefail; \
    bcftools concat -f '${TMP}/concat_list.txt' --threads ${T} \
        -Oz -o '${OUT}/allsites.hc.vcf.gz'"
run tabix -p vcf "${OUT}/allsites.hc.vcf.gz"

N=$(run bash -c "bcftools index -n '${OUT}/allsites.hc.vcf.gz'")
echo "[$(date)] all-sites sites : ${N}"
