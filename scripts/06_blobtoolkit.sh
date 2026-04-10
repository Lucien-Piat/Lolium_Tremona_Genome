#!/bin/bash
#SBATCH --job-name=blobtoolkit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/06_blobtoolkit_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/blobtoolkit.sif)
READS="raw_reads/lmultiflorum_hifi.fastq.gz"
ASM_GZ="results/04_organellar/lmultiflorum.nuclear.fa.gz"
OUTDIR="results/05_blobtoolkit"
TRACKING="results/assembly_tracking.tsv"
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)

TAXDUMP="reference_data/taxdump"
BLASTDB_DIR="/cluster/project/clcgenomics/CLC_BLAST_DB"
BLASTDB_NAME="nt"
TAXID=4521

run() { singularity exec --bind "${BLASTDB_DIR}" --env BLASTDB="${BLASTDB_DIR}" "${SIF}" "$@"; }

track() {
    local stage=$1 fa=$2
    local stats nseq size
    stats=$(run seqkit stats -T "${fa}" | tail -1)
    nseq=$(echo "${stats}" | cut -f4)
    size=$(echo "${stats}" | cut -f5)
    printf '%s\t%s\t%s\t%s\n' "${stage}" "$(readlink -f "${fa}")" "${nseq}" "${size}" >> "${TRACKING}"
    echo "${stage}: ${nseq} contigs, ${size} bp"
}

mkdir -p "${OUTDIR}" logs
[[ -f "${TRACKING}" ]] || printf 'stage\tfile\tcontigs\tsize\n' > "${TRACKING}"

run pigz -dcp "${T}" "${ASM_GZ}" > "${OUTDIR}/assembly.fa"

echo "Mapping reads..."
run minimap2 -t "${T}" -ax map-hifi "${OUTDIR}/assembly.fa" "${READS}" \
    | run samtools sort -@ "${T}" -O BAM -o "${OUTDIR}/coverage.bam" -
run samtools index -@ "${T}" "${OUTDIR}/coverage.bam"

echo "Running blastn against nt..."
run blastn \
    -db "${BLASTDB_DIR}/${BLASTDB_NAME}" \
    -query "${OUTDIR}/assembly.fa" \
    -outfmt "6 qseqid staxids bitscore std" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -evalue 1e-25 \
    -num_threads "${T}" \
    -out "${OUTDIR}/blast.out"

echo "Creating BlobDir..."
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

echo "Filtering contigs..."
run blobtools filter \
    --param bestsumorder_phylum--Keys=Streptophyta \
    --param bestsumorder_phylum--Keys=no-hit \
    --fasta "${OUTDIR}/assembly.fa" \
    --output "${OUTDIR}/filtered" \
    "${OUTDIR}/blobdir"

FILTERED="${OUTDIR}/filtered/assembly.filtered.fa"
if [[ -f "${FILTERED}" ]]; then
    mv "${FILTERED}" "${OUTDIR}/lmultiflorum.decontam.fa"
else
    echo "WARN: no filtered output, keeping assembly as-is"
    cp "${OUTDIR}/assembly.fa" "${OUTDIR}/lmultiflorum.decontam.fa"
fi

run pigz -p "${T}" "${OUTDIR}/lmultiflorum.decontam.fa"

# ---- Cleanup ----
rm -f "${OUTDIR}/assembly.fa" "${OUTDIR}/coverage.bam" "${OUTDIR}/coverage.bam.bai"
rm -f "${OUTDIR}/blast.out"
rm -rf "${OUTDIR}/filtered"

cd "${ROOT}"
track "assembly-decontam" "${OUTDIR}/lmultiflorum.decontam.fa.gz"

echo "BlobDir kept at: ${OUTDIR}/blobdir (use 'blobtools view' to explore)"