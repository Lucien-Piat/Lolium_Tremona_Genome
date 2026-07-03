#!/bin/bash
#SBATCH --job-name=jc_prep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/04_prep_%j.log

set -euo pipefail
cd "$(readlink -f .)"
N=${N:-100}
SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
GVCFDIR=$(readlink -f results/03_gvcf)
OUTDIR="results/04_joint_calling"
MAP="${OUTDIR}/sample_map.tsv"
mkdir -p "${OUTDIR}" logs
run() { singularity exec --bind "$PWD" "${SIF}" "$@"; }

REF=$(readlink -f "$(run yq -r '.reference_genome' "${YAML}")")
[[ -s "${REF}" && -s "${REF}.fai" && -s "${REF%.fa}.dict" ]] \
    || { echo "ERROR: reference, .fai ou .dict manquant" >&2; exit 1; }

for g in "${GVCFDIR}"/*.g.vcf.gz; do
    s=$(basename "${g}" .g.vcf.gz)
    printf "%s\t%s\n" "${s}" "$(readlink -f "${g}")"
done > "${MAP}"

miss=0
while read -r _ g; do [[ -s "${g}.tbi" ]] || { echo "index manquant: ${g}" >&2; miss=1; }; done < "${MAP}"
[[ "${miss}" -eq 0 ]] || { echo "ERROR: des .tbi manquent" >&2; exit 1; }

rm -rf "${OUTDIR}/scatter"
run gatk SplitIntervals -R "${REF}" --scatter-count "${N}" -O "${OUTDIR}/scatter"
