#!/bin/bash
#SBATCH --job-name=blobtoolkit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/06_blobtoolkit_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/blobtoolkit.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
ASM_GZ=$(readlink -f results/04_organellar/lmultiflorum.nuclear.fa.gz)
OUTDIR="results/05_blobtoolkit"
TRACKING="results/assembly_tracking.tsv"
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)

TAXDUMP=$(readlink -f reference_data/ncbi_taxdump)
BLASTDB_DIR="/cluster/project/clcgenomics/CLC_BLAST_DB"
BLASTDB_NAME="nt"
TAXID=4521

# Bind both scratch (for project files) and the cluster project dir (for BLASTDB).
BIND="/cluster/scratch,${BLASTDB_DIR}"

CACHE="$(pwd)/.cache_blobtoolkit"
mkdir -p "${CACHE}/home" "${CACHE}/matplotlib" "${CACHE}/fontconfig" \
         "${OUTDIR}" "${TAXDUMP}" logs

run() {
    singularity exec \
        --bind "${BIND}" \
        --env BLASTDB="${BLASTDB_DIR}" \
        --env HOME="${CACHE}/home" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env FONTCONFIG_PATH="${CACHE}/fontconfig" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" "$@"
}

run_sh() {
    singularity exec \
        --bind "${BIND}" \
        --env BLASTDB="${BLASTDB_DIR}" \
        --env HOME="${CACHE}/home" \
        --env MPLCONFIGDIR="${CACHE}/matplotlib" \
        --env FONTCONFIG_PATH="${CACHE}/fontconfig" \
        --env XDG_CACHE_HOME="${CACHE}" \
        "${SIF}" bash -c "$@"
}

track() {
    local stage=$1 fa=$2
    local stats nseq size
    stats=$(run seqkit stats -T "${fa}" | tail -1)
    nseq=$(echo "${stats}" | cut -f4)
    size=$(echo "${stats}" | cut -f5)
    printf '%s\t%s\t%s\t%s\n' "${stage}" "$(readlink -f "${fa}")" "${nseq}" "${size}" >> "${TRACKING}"
    echo "${stage}: ${nseq} contigs, ${size} bp"
}

[[ -f "${TRACKING}" ]] || printf 'stage\tfile\tcontigs\tsize\n' > "${TRACKING}"

if [[ ! -f "${TAXDUMP}/nodes.dmp" ]]; then
    echo "Downloading NCBI taxdump..."
    wget -q -O "${TAXDUMP}/new_taxdump.tar.gz" \
        https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
    tar xzf "${TAXDUMP}/new_taxdump.tar.gz" -C "${TAXDUMP}"
    rm -f "${TAXDUMP}/new_taxdump.tar.gz"
else
    echo "Using cached taxdump at ${TAXDUMP}"
fi

run pigz -dcp 4 "${ASM_GZ}" > "${OUTDIR}/assembly.fa"

echo "Step 1: mapping HiFi reads with minimap2"
run_sh "minimap2 -t $((T-4)) -ax map-hifi ${OUTDIR}/assembly.fa ${READS} \
         | samtools sort -@ 4 -m 2G -O BAM -o ${OUTDIR}/coverage.bam -"
run samtools index -@ 4 "${OUTDIR}/coverage.bam"

echo "Step 2: blastn against nt (this is the long step)"
run blastn \
    -db "${BLASTDB_DIR}/${BLASTDB_NAME}" \
    -query "${OUTDIR}/assembly.fa" \
    -outfmt "6 qseqid staxids bitscore std" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -evalue 1e-25 \
    -num_threads "${T}" \
    -out "${OUTDIR}/blast.out"

echo "Step 3: building BlobDir"
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

echo "Step 4: filtering to keep Streptophyta + no-hit contigs"
run blobtools filter \
    --param bestsumorder_phylum--Keys=Streptophyta \
    --param bestsumorder_phylum--Keys=no-hit \
    --fasta "${OUTDIR}/assembly.fa" \
    --output "${OUTDIR}/filtered" \
    "${OUTDIR}/blobdir"

FILTERED="${OUTDIR}/filtered/assembly.filtered.fa"
if [[ -f "${FILTERED}" ]]; then
    mv "${FILTERED}" "${OUTDIR}/lmultiflorum.decontam.fa"
    echo "INFO: filtered assembly produced"
else
    echo "WARN: no filtered output, keeping assembly as-is"
    mv "${OUTDIR}/assembly.fa" "${OUTDIR}/lmultiflorum.decontam.fa"
fi

if run pigz -p "${T}" "${OUTDIR}/lmultiflorum.decontam.fa"; then
    rm -f "${OUTDIR}/assembly.fa" \
          "${OUTDIR}/coverage.bam" "${OUTDIR}/coverage.bam.bai" \
          "${OUTDIR}/blast.out"
    rm -rf "${OUTDIR}/filtered" "${CACHE}"
else
    echo "ERROR: output compression failed, keeping intermediates for debug" >&2
    exit 1
fi

cd "${ROOT}"
track "assembly-decontam" "${OUTDIR}/lmultiflorum.decontam.fa.gz"