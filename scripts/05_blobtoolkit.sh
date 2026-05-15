#!/bin/bash
#SBATCH --job-name=blobtoolkit
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=36:00:00
#SBATCH --output=logs/05_blobtoolkit_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/blobtoolkit.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
ASM=$(readlink -f results/04c_purgegrass/final_primary_with_trim.fa)
TAXDUMP=$(readlink -f reference_data/ncbi_taxdump)

OUTDIR="results/05_blobtoolkit_purged"
FINAL_BASENAME="lmultiflorum.purgegrass.decontam.fa"
BLASTDB_DIR="/cluster/project/clcgenomics/CLC_BLAST_DB"
TAXID=4521

T=${SLURM_CPUS_PER_TASK:-4}
BIND="/cluster/scratch,${BLASTDB_DIR}"
CACHE="$(pwd)/.cache_btk"

mkdir -p "${CACHE}" "${OUTDIR}" logs

# Setup singularity environment natively
export SINGULARITYENV_BLASTDB="${BLASTDB_DIR}"
export SINGULARITYENV_MPLCONFIGDIR="${CACHE}"
export SINGULARITYENV_FONTCONFIG_PATH="${CACHE}"
export SINGULARITYENV_XDG_CACHE_HOME="${CACHE}"

run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"; }

# Minimal safety checks
for f in "${SIF}" "${READS}" "${ASM}"; do
    [[ -s "$f" ]] || { echo "ERROR: Missing $f" >&2; exit 1; }
done

# Stage assembly
ASM_LOCAL="${OUTDIR}/assembly.fa"
ln -sf "${ASM}" "${ASM_LOCAL}"
[[ -s "${ASM_LOCAL}.fai" ]] || run samtools faidx "${ASM_LOCAL}"

# Step 1: Mapping (resume if BAM index exists)
COVBAM="${OUTDIR}/coverage.bam"
[[ -s "${COVBAM}" ]] && [[ ! -s "${COVBAM}.bai" ]] && rm -f "${COVBAM}" # clean stale bam

if [[ ! -s "${COVBAM}.bai" ]]; then
    run_sh "minimap2 -t $((T-2)) -I 4G -K 5G -ax map-hifi --secondary=no ${ASM_LOCAL} ${READS} \
            | samtools sort -@ 2 -m 2G -O BAM -o ${COVBAM} -"
    run samtools index -@ "${T}" "${COVBAM}"
fi

# Step 2: BLAST (resume if output exists)
BLASTOUT="${OUTDIR}/blast.out"
if [[ ! -s "${BLASTOUT}" ]]; then
    run blastn -db nt -query "${ASM_LOCAL}" -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads "${T}" -out "${BLASTOUT}"
fi

# Step 3: Build BlobDir (resume if meta.json exists)
BLOBDIR="${OUTDIR}/blobdir"
if [[ ! -f "${BLOBDIR}/meta.json" ]]; then
    run blobtools create --fasta "${ASM_LOCAL}" --taxid "${TAXID}" --taxdump "${TAXDUMP}" "${BLOBDIR}"
    run blobtools add --cov "${COVBAM}" "${BLOBDIR}"
    run blobtools add --hits "${BLASTOUT}" --taxrule bestsumorder --taxdump "${TAXDUMP}" "${BLOBDIR}"
fi

# Step 4: Filter (Streptophyta + no-hit)
rm -rf "${OUTDIR}/filtered"
run blobtools filter \
    --param bestsumorder_phylum--Keys=Streptophyta \
    --param bestsumorder_phylum--Keys=no-hit \
    --fasta "${ASM_LOCAL}" \
    --output "${OUTDIR}/filtered" \
    "${BLOBDIR}"

# Finalize and compress
FINAL="${OUTDIR}/${FINAL_BASENAME}"
FILTERED=$(find "${OUTDIR}/filtered" \( -name "*.filtered.fasta" -o -name "*.filtered.fa" \) 2>/dev/null | head -n 1)

if [[ -n "${FILTERED}" ]] && [[ -s "${FILTERED}" ]]; then
    mv "${FILTERED}" "${FINAL}"
else
    cp "${ASM_LOCAL}" "${FINAL}"
fi

run pigz -f -p "${T}" "${FINAL}"

# Cleanup
rm -f "${ASM_LOCAL}" "${ASM_LOCAL}.fai" "${COVBAM}" "${COVBAM}.bai" "${BLASTOUT}"
rm -rf "${OUTDIR}/filtered" "${CACHE}"