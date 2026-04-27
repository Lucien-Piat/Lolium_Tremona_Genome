#!/bin/bash
#SBATCH --job-name=purge_dups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=2G
#SBATCH --time=15:00:00
#SBATCH --output=logs/04_purge_dups_%j.log

set -euo pipefail

SIF=$(readlink -f images/sif/polish.sif)
READS=$(readlink -f raw_reads/lmultiflorum_hifi.fastq.gz)
ASM_GZ=$(readlink -f results/02_assembly/lmultiflorum.bp.p_ctg.fa.gz)
OUTDIR="results/03_purge"
TRACKING="results/assembly_tracking.tsv"
T=${SLURM_CPUS_PER_TASK:-4}
ROOT=$(pwd)

BIND="/cluster/scratch"

run()    { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"; }

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

run pigz -dcp 4 "${ASM_GZ}" > "${OUTDIR}/assembly.fa"
cd "${OUTDIR}"

echo "Step 1: reads vs assembly (minimap2 map-hifi)"
run_sh "minimap2 -t $((T-4)) -xmap-hifi assembly.fa ${READS} | pigz -p 4 > aligned.paf.gz"

echo "Step 2: coverage stats and cutoffs"
run pbcstat aligned.paf.gz
run calcuts PB.stat > cutoffs

echo "Step 3: self-alignment"
run split_fa assembly.fa > asm.split.fa
run_sh "minimap2 -t $((T-4)) -xasm5 -DP asm.split.fa asm.split.fa | pigz -p 4 > self_aln.paf.gz"

echo "Step 4: purge"
run purge_dups -2 -T cutoffs -c PB.base.cov assembly.fa > dups.bed
run get_seqs -e dups.bed assembly.fa || true

if [[ -s purged.fa ]]; then
    mv purged.fa lmultiflorum.purged.fa
    [[ -f hap.fa ]] && mv hap.fa lmultiflorum.haplotigs.fa || touch lmultiflorum.haplotigs.fa
else
    echo "WARN: nothing purged, keeping original assembly"
    mv assembly.fa lmultiflorum.purged.fa
    touch lmultiflorum.haplotigs.fa
fi

if run pigz -p "${T}" lmultiflorum.purged.fa lmultiflorum.haplotigs.fa; then
    rm -f aligned.paf.gz self_aln.paf.gz asm.split.fa assembly.fa \
          PB.stat PB.base.cov PB.cov2.bin PB.cov.wig cutoffs dups.bed
else
    echo "ERROR: output compression failed, keeping intermediates for debug" >&2
    exit 1
fi

cd "${ROOT}"
track "assembly-purged" "${OUTDIR}/lmultiflorum.purged.fa.gz"