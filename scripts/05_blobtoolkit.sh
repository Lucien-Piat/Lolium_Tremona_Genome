#!/bin/bash
#SBATCH --job-name=blobtoolkit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=36:00:00
#SBATCH --output=logs/05_blobtoolkit_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/blobtoolkit.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
ASM_GZ=$(readlink -f results/04c_purgegrass/lmultiflorum.purgegrass.fa.gz)
OUTDIR="results/05_blobtoolkit_purged"
T=${SLURM_CPUS_PER_TASK:-4}

TAXDUMP=$(readlink -f reference_data/ncbi_taxdump)
BLASTDB_DIR="/cluster/project/clcgenomics/CLC_BLAST_DB"
BLASTDB_NAME="nt"
TAXID=4521

BIND="/cluster/scratch,${BLASTDB_DIR}"

CACHE="$(pwd)/.cache_btk_purged"
mkdir -p "${CACHE}/home" "${CACHE}/matplotlib" "${CACHE}/fontconfig" "${OUTDIR}" logs

run() {
    singularity exec \
        --bind "${BIND}" \
        --env BLASTDB="${BLASTDB_DIR}" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env FONTCONFIG_PATH="${CACHE}/fontconfig" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

run_sh() {
    singularity exec \
        --bind "${BIND}" \
        --env BLASTDB="${BLASTDB_DIR}" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env FONTCONFIG_PATH="${CACHE}/fontconfig" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" bash -c "$@"
}

[[ -f "${TAXDUMP}/nodes.dmp" ]] || \
    { echo "ERROR: taxdump missing at ${TAXDUMP}" >&2; exit 1; }

echo "==================================================="
echo "BlobToolKit on purged assembly"
echo "Started: $(date)"
echo "Resources: ${T} cpus, $((T*4)) GB RAM"
echo "==================================================="

# Decompress assembly
run pigz -dcp "${T}" "${ASM_GZ}" > "${OUTDIR}/assembly.fa"

# Step 1: map HiFi reads to assembly
# -I 2G keeps memory predictable within cgroup limits
echo "Step 1: minimap2 mapping at $(date)"
run_sh "minimap2 -t $((T-2)) -I 2G -K 5G -ax map-hifi --secondary=no \
          ${OUTDIR}/assembly.fa ${READS} \
        | samtools sort -@ 2 -m 2G -O BAM -o ${OUTDIR}/coverage.bam -"
run samtools index -@ "${T}" "${OUTDIR}/coverage.bam"

# Step 2: blastn against nt (this is the long step, hours)
echo "Step 2: blastn at $(date)"
run blastn \
    -db "${BLASTDB_DIR}/${BLASTDB_NAME}" \
    -query "${OUTDIR}/assembly.fa" \
    -outfmt "6 qseqid staxids bitscore std" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -evalue 1e-25 \
    -num_threads "${T}" \
    -out "${OUTDIR}/blast.out"

# Step 3: build BlobDir
echo "Step 3: building BlobDir at $(date)"
run blobtools create \
    --fasta "${OUTDIR}/assembly.fa" \
    --taxid "${TAXID}" \
    --taxdump "${TAXDUMP}" \
    "${OUTDIR}/blobdir"

run blobtools add \
    --cov "${OUTDIR}/coverage.bam" \
    "${OUTDIR}/blobdir"

run blobtools add \
    --hits "${OUTDIR}/blast.out" \
    --taxrule bestsumorder \
    --taxdump "${TAXDUMP}" \
    "${OUTDIR}/blobdir"

# Step 4: filter to keep Streptophyta and no-hit contigs
echo "Step 4: filtering at $(date)"
run blobtools filter \
    --param bestsumorder_phylum--Keys=Streptophyta \
    --param bestsumorder_phylum--Keys=no-hit \
    --fasta "${OUTDIR}/assembly.fa" \
    --output "${OUTDIR}/filtered" \
    "${OUTDIR}/blobdir"

FILTERED="${OUTDIR}/filtered/assembly.filtered.fa"
if [[ -f "${FILTERED}" ]]; then
    mv "${FILTERED}" "${OUTDIR}/lmultiflorum.purged.decontam.fa"
    echo "INFO: filtered assembly produced"
else
    echo "WARN: no filtered output, keeping assembly as-is"
    cp "${OUTDIR}/assembly.fa" "${OUTDIR}/lmultiflorum.purged.decontam.fa"
fi

# Compress and clean up heavy intermediates
run pigz -p "${T}" "${OUTDIR}/lmultiflorum.purged.decontam.fa"

rm -f "${OUTDIR}/assembly.fa" \
      "${OUTDIR}/coverage.bam" "${OUTDIR}/coverage.bam.bai" \
      "${OUTDIR}/blast.out"
rm -rf "${OUTDIR}/filtered" "${CACHE}"

echo ""
echo "==================================================="
echo "BlobToolKit done at $(date)"
echo "==================================================="
echo "Filtered assembly: ${OUTDIR}/lmultiflorum.purged.decontam.fa.gz"
echo "BlobDir for plots: ${OUTDIR}/blobdir/"