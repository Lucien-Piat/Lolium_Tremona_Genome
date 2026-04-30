#!/bin/bash
#SBATCH --job-name=pg_prep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/04c_pg_01_prep_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/purgegrass.sif)
ASM_GZ=$(readlink -f results/02_assembly/S_lmultiflorum.bp.p_ctg.fa.gz)
REF_FA=$(readlink -f reference_data/full_ref_genome/LOLMU.fa)
REF_GFF=$(readlink -f reference_data/full_ref_genome/LOLMU.gff)
OUTDIR="results/04c_purgegrass"
T=${SLURM_CPUS_PER_TASK:-4}

BIND=$PWD
run() { singularity exec --bind "${BIND}" "${SIF}" "$@"; }

mkdir -p "${OUTDIR}" reference_data/transcripts logs

# Extract transcripts from the reference annotation
TRANSCRIPTS="reference_data/transcripts/lolmu_transcripts.fa"
if [[ ! -f "${TRANSCRIPTS}" ]]; then
    echo "Extracting transcripts from ${REF_GFF}..."
    run gffread -w "${TRANSCRIPTS}" -g "${REF_FA}" "${REF_GFF}"
    run seqkit stats "${TRANSCRIPTS}"
fi

# Decompress assembly and create .fai index
ASM_FA="${OUTDIR}/assembly.fa"
if [[ ! -f "${ASM_FA}" ]]; then
    echo "Decompressing assembly..."
    run pigz -dcp "${T}" "${ASM_GZ}" > "${ASM_FA}"
    run samtools faidx "${ASM_FA}"
fi

echo "Prep done at $(date)"
ls -lh "${OUTDIR}/"
ls -lh reference_data/transcripts/