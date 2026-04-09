#!/bin/bash
#SBATCH --job-name=purge_dups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/04_purge_dups_%j.log

SIF="$(pwd)/images/sif/polish.sif"
READS="$(pwd)/raw_reads/lmultiflorum_hifi.fastq.gz"
ASM_GZ="$(pwd)/results/02_assembly/lmultiflorum.bp.p_ctg.fa.gz"
OUTDIR="$(pwd)/results/03_purge"
THREADS=${SLURM_CPUS_PER_TASK}

mkdir -p "${OUTDIR}" logs

singularity exec "${SIF}" pigz -dk -p "${THREADS}" "${ASM_GZ}"
ASM="${ASM_GZ%.gz}"

cd "${OUTDIR}"

# Map HiFi reads to assembly
singularity exec "${SIF}" \
    minimap2 -t "${THREADS}" -xmap-hifi "${ASM}" "${READS}" \
    | singularity exec "${SIF}" pigz -p 4 \
    > aligned.paf.gz

# Coverage stats and cutoffs
singularity exec "${SIF}" pbcstat aligned.paf.gz
singularity exec "${SIF}" calcuts PB.stat > cutoffs

# Self-alignment
singularity exec "${SIF}" split_fa "${ASM}" > asm.split.fa
singularity exec "${SIF}" \
    minimap2 -t "${THREADS}" -xasm5 -DP asm.split.fa asm.split.fa \
    | singularity exec "${SIF}" pigz -p 4 \
    > self_aln.paf.gz

# Purge
singularity exec "${SIF}" \
    purge_dups -2 -T cutoffs -c PB.base.cov "${ASM}" > dups.bed

# Get purged + haplotig sequences
singularity exec "${SIF}" get_seqs -e dups.bed "${ASM}"
mv purged.fa lmultiflorum.purged.fa
mv hap.fa lmultiflorum.haplotigs.fa

echo "=== PRE-PURGE ==="
singularity exec "${SIF}" seqkit stats -a "${ASM}"
echo "=== POST-PURGE ==="
singularity exec "${SIF}" seqkit stats -a lmultiflorum.purged.fa
echo "=== HAPLOTIGS ==="
singularity exec "${SIF}" seqkit stats -a lmultiflorum.haplotigs.fa

singularity exec "${SIF}" pigz -p "${THREADS}" lmultiflorum.purged.fa
singularity exec "${SIF}" pigz -p "${THREADS}" lmultiflorum.haplotigs.fa

rm -f aligned.paf.gz self_aln.paf.gz asm.split.fa
rm -f PB.stat PB.base.cov PB.cov2.bin cutoffs dups.bed
rm -f "${ASM}"