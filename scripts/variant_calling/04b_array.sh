#!/bin/bash
#SBATCH --job-name=jc_shard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=20:00:00
#SBATCH --array=0-199%25
#SBATCH --output=logs/04_shard_%A_%a.log

# Step 2: import + genotype one shard per array task.
# Run AFTER prep: sbatch scripts/04b_shards.sh

set -euo pipefail
cd "$(readlink -f .)"
SIF=$(readlink -f images/sif/varcall.sif)
YAML=$(readlink -f reads/samples.yaml)
OUTDIR="results/04_joint_calling"
MAP="${OUTDIR}/sample_map.tsv"
T=${SLURM_CPUS_PER_TASK:-1}
run() { singularity exec --bind "$PWD" "${SIF}" "$@"; }

[[ -s "${MAP}" ]] || { echo "ERROR: ${MAP} manquant, lance 04a_prep.sh d'abord" >&2; exit 1; }
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