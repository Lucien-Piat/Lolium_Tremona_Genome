#!/bin/bash
#SBATCH --job-name=ld_decay
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/ld_%x_%j.log
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/popgen.sif"
VCF="${ROOT}/results/filtered_vcf/snps.masked.biallelic.vcf.gz"
POP="${ROOT}/scripts/pop_gen/pop.tsv"
OUT="${ROOT}/results/ld_decay"

POPS="GULF L31 L46 L60 PR SLB TREM CH"

export APPTAINER_HOME="${ROOT}/.apphome"
export SINGULARITY_HOME="${APPTAINER_HOME}"
mkdir -p "${APPTAINER_HOME}" "${OUT}"

run() { singularity exec --home "${APPTAINER_HOME}" --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

NMIN=9999
for p in ${POPS}; do
    n=$(awk -v pp="${p}" '$2==pp' "${POP}" | wc -l)
    echo "[$(date)] ${p} : ${n} individus"
    [ "${n}" -lt "${NMIN}" ] && NMIN=${n}
done

# LD decay par population
for p in ${POPS}; do
    run bash -c "awk -v pp='${p}' '\$2==pp{print \$1}' '${POP}' \
        | shuf | head -n ${NMIN} > '${OUT}/${p}.txt'"
    echo "[$(date)] LD decay : ${p} (n=${NMIN})"
    run PopLDdecay -InVCF "${VCF}" \
        -SubPop "${OUT}/${p}.txt" \
        -MaxDist 300 \
        -OutStat "${OUT}/${p}.stat"
done

run bash -c "for p in ${POPS}; do echo -e '${OUT}/'\$p'.stat.gz\t'\$p; done > '${OUT}/multi.list'"
