#!/bin/bash
# Controller: run from the login node with `bash scripts/04_joint_call.sh`
# Chains 3 SLURM jobs: prep -> shard array -> gather.
set -euo pipefail

PROJECT=$(readlink -f .)
N=100                      # number of shards (also the array size)
CONC=40                    # max shards running at once
OUTDIR="${PROJECT}/results/04_joint_calling"
mkdir -p "${OUTDIR}" "${PROJECT}/logs"

COMMON_EXPORT="ALL,PROJECT=${PROJECT},N=${N}"

PREP_ID=$(sbatch --parsable \
    --job-name=jc_prep \
    --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=00:20:00 \
    --output="${PROJECT}/logs/04_prep_%j.log" \
    --export="${COMMON_EXPORT}" <<'EOF'
#!/bin/bash
set -euo pipefail
cd "${PROJECT}"
SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
GVCFDIR=$(readlink -f results/03_gvcf)
OUTDIR="results/04_joint_calling"
MAP="${OUTDIR}/sample_map.tsv"
run() { singularity exec --bind "$PWD" "${SIF}" "$@"; }

REF=$(readlink -f "$(run yq -r '.reference_genome' "${YAML}")")
[[ -s "${REF}" && -s "${REF}.fai" && -s "${REF%.fa}.dict" ]] \
    || { echo "ERROR: reference, .fai ou .dict manquant" >&2; exit 1; }

for g in "${GVCFDIR}"/*.g.vcf.gz; do
    s=$(basename "${g}" .g.vcf.gz)
    printf "%s\t%s\n" "${s}" "$(readlink -f "${g}")"
done > "${MAP}"
echo "Cohorte : $(wc -l < "${MAP}") echantillons"

miss=0
while read -r _ g; do [[ -s "${g}.tbi" ]] || { echo "index manquant: ${g}" >&2; miss=1; }; done < "${MAP}"
[[ "${miss}" -eq 0 ]] || { echo "ERROR: des .tbi manquent" >&2; exit 1; }

rm -rf "${OUTDIR}/scatter"
run gatk SplitIntervals -R "${REF}" --scatter-count "${N}" -O "${OUTDIR}/scatter"
echo "Shards crees : $(ls "${OUTDIR}/scatter"/*-scattered.interval_list | wc -l)"
EOF
)
echo "prep   submitted : ${PREP_ID}"

# 2. SHARDS
ARR_ID=$(sbatch --parsable \
    --job-name=jc_shard \
    --dependency=afterok:"${PREP_ID}" \
    --array=0-$((N-1))%"${CONC}" \
    --ntasks=1 --cpus-per-task=1 --mem-per-cpu=20G --time=06:00:00 \
    --output="${PROJECT}/logs/04_shard_%A_%a.log" \
    --export="${COMMON_EXPORT}" <<'EOF'
#!/bin/bash
set -euo pipefail
cd "${PROJECT}"
SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
OUTDIR="results/04_joint_calling"
MAP="${OUTDIR}/sample_map.tsv"
T=${SLURM_CPUS_PER_TASK:-1}
run() { singularity exec --bind "$PWD" "${SIF}" "$@"; }

REF=$(readlink -f "$(run yq -r '.reference_genome' "${YAML}")")
ID=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
IL="${OUTDIR}/scatter/${ID}-scattered.interval_list"
[[ -s "${IL}" ]] || { echo "ERROR: ${IL} manquant" >&2; exit 1; }

SHARD_OUT="${OUTDIR}/shards/${ID}"
WORKTMP="${SHARD_OUT}/tmp"
DB="${WORKTMP}/genomicsdb"
mkdir -p "${SHARD_OUT}"
rm -rf "${WORKTMP}"; mkdir -p "${WORKTMP}"
trap 'rm -rf "${WORKTMP}"' EXIT

# heaps sized to a 22.5Mb region: 8G import, 12G genotype, ~8G left for
# GenomicsDB native memory inside the 20G allocation.
run gatk --java-options "-Xmx8G" GenomicsDBImport \
    --sample-name-map "${MAP}" \
    --genomicsdb-workspace-path "${DB}" \
    -L "${IL}" \
    --merge-input-intervals \
    --batch-size 50 \
    --reader-threads "${T}" \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir "${WORKTMP}"

run gatk --java-options "-Xmx12G" GenotypeGVCFs \
    -R "${REF}" \
    -V "gendb://${DB}" \
    -L "${IL}" \
    --include-non-variant-sites \
    -O "${SHARD_OUT}/allsites.vcf.gz" \
    --tmp-dir "${WORKTMP}"

[[ -s "${SHARD_OUT}/allsites.vcf.gz" && -s "${SHARD_OUT}/allsites.vcf.gz.tbi" ]] \
    || { echo "ERROR: shard ${ID} incomplet" >&2; exit 1; }
echo "shard ${ID} OK : $(run bcftools index -n "${SHARD_OUT}/allsites.vcf.gz") sites"
EOF
)
echo "shards submitted : ${ARR_ID}"


GAT_ID=$(sbatch --parsable \
    --job-name=jc_gather \
    --dependency=afterok:"${ARR_ID}" \
    --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8G --time=06:00:00 \
    --output="${PROJECT}/logs/04_gather_%j.log" \
    --export="${COMMON_EXPORT}" <<'EOF'
#!/bin/bash
set -euo pipefail
cd "${PROJECT}"
SIF=$(readlink -f images/sif/varcall.sif)
OUTDIR="results/04_joint_calling"
run()    { singularity exec --bind "$PWD" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "$PWD" "${SIF}" bash -c "$@"; }

got=$(ls "${OUTDIR}"/shards/*/allsites.vcf.gz 2>/dev/null | wc -l)
[[ "${got}" -eq "${N}" ]] || { echo "ERROR: ${got}/${N} shards seulement" >&2; exit 1; }

ls "${OUTDIR}"/shards/*/allsites.vcf.gz | sort > "${OUTDIR}/vcf.list"
rm -f "${OUTDIR}/cohort_allsites.vcf.gz" "${OUTDIR}/cohort_allsites.vcf.gz.tbi"

ALLSITES="${OUTDIR}/cohort_allsites.vcf.gz"
# block-copy shards (no recompression); fall back to threaded recompress if needed
if ! run_sh "bcftools concat --naive-force -f '${OUTDIR}/vcf.list' -o '${ALLSITES}'"; then
    echo "naive-force indisponible, recompression (--threads 2)" >&2
    run_sh "bcftools concat --threads 2 -f '${OUTDIR}/vcf.list' -Oz -o '${ALLSITES}'"
fi
run bcftools index -t "${ALLSITES}"
N_ALL=$(run bcftools index -n "${ALLSITES}")
echo "all-sites       : ${ALLSITES} (${N_ALL} sites)"

# variants-only : SNPs, QUAL>20, QD>8 (Stritt et al. 2022)
SNPS="${OUTDIR}/cohort_snps_filtered.vcf.gz"
run_sh "bcftools view -v snps '${ALLSITES}' \
        | bcftools filter -i 'QUAL>20 && INFO/QD>8' --threads 2 -Oz -o '${SNPS}'"
run bcftools index -t "${SNPS}"
N_SNP=$(run bcftools index -n "${SNPS}")
run bcftools stats "${SNPS}" > "${SNPS%.vcf.gz}.stats"
echo "variants filtre : ${SNPS} (${N_SNP} SNPs)"

# cleanup: outputs verified by index, free ~200G of redundant shards
if [[ "${N_ALL}" -gt 0 && "${N_SNP}" -gt 0 && -s "${ALLSITES}.tbi" && -s "${SNPS}.tbi" ]]; then
    du -sh "${OUTDIR}/shards" 2>/dev/null | awk '{print "shards liberes : "$1}'
    rm -rf "${OUTDIR}/shards"
    rm -f "${OUTDIR}/vcf.list"
    echo "shards supprimes (verif OK)"
else
    echo "WARN: verification incomplete, shards CONSERVES" >&2
    exit 1
fi
EOF
)
echo "gather submitted : ${GAT_ID}"

echo
echo "Chaine : ${PREP_ID} -> ${ARR_ID} (array 0-$((N-1))%${CONC}) -> ${GAT_ID}"