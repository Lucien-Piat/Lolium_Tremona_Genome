#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3G
#SBATCH --time=24:00:00
#SBATCH --output=logs/03_hifiasm_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/hifiasm.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
OUTDIR="results/02_assembly"
TRACKING="results/assembly_tracking.tsv"
T=${SLURM_CPUS_PER_TASK}
ROOT=$(pwd)

run() { singularity exec "${SIF}" "$@"; }

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

ln -sf "${READS}" "${OUTDIR}/reads.fastq.gz"
cd "${OUTDIR}"

run hifiasm -t "${T}" -o lmultiflorum --telo-m TTTAGGG reads.fastq.gz

for TYPE in bp.p_ctg bp.hap1.p_ctg bp.hap2.p_ctg; do
    awk '/^S/{print ">"$2; print $3}' "lmultiflorum.${TYPE}.gfa" \
        | run pigz -p "${T}" > "lmultiflorum.${TYPE}.fa.gz"
done

rm -f lmultiflorum*.bin reads.fastq.gz
run pigz -p "${T}" lmultiflorum*.gfa

cd "${ROOT}"
track "assembly-raw"  "${OUTDIR}/lmultiflorum.bp.p_ctg.fa.gz"
track "hap1-raw"      "${OUTDIR}/lmultiflorum.bp.hap1.p_ctg.fa.gz"
track "hap2-raw"      "${OUTDIR}/lmultiflorum.bp.hap2.p_ctg.fa.gz"