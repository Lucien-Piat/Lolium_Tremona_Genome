#!/bin/bash
# Build distance tree from BCALM2 assemblies + reference genomes.
#SBATCH --job-name=mashtree
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=logs/mashtree_%j.log

set -euo pipefail
exec 2>&1
SIF="${SIF_MASHTREE:-./images/sif/draft_mash_tree.sif}"
RESULTS="${RESULTS:-./results}"
ASSEMBLIES_DIR="${ASSEMBLIES_DIR:-./assemblies}"
GENOMES_DIR="${GENOMES_DIR:-$(cd ../genomes 2>/dev/null && pwd || echo '')}"
export TMPDIR="$(pwd)/tmp_mash"
mkdir -p "$RESULTS" "$TMPDIR"

N_ASM=$(ls -1 "${ASSEMBLIES_DIR}"/*.unitigs.fa.gz 2>/dev/null | wc -l)
N_GEN=$(ls -1 "${GENOMES_DIR}"/*.fasta.gz 2>/dev/null | wc -l)
echo "[INFO] $N_ASM assemblies, $N_GEN reference genomes"

apptainer exec \
    --no-home \
    --pwd "$(pwd)" \
    --bind "$(pwd):$(pwd),${GENOMES_DIR}:${GENOMES_DIR}" \
    "$SIF" \
    mashtree \
        --numcpus 8 \
        --outtree "${RESULTS}/mashtree.dnd" \
        --outmatrix "${RESULTS}/distances.tsv" \
        --kmerlength 21 \
        --sketch-size 10000 \
        --mindepth 0 \
        "${ASSEMBLIES_DIR}"/*.unitigs.fa.gz \
        "${GENOMES_DIR}"/*.fasta.gz

