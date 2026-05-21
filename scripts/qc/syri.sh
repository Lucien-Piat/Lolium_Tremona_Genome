#!/bin/bash
#SBATCH --job-name=syri_tremona
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --time=12:00:00
#SBATCH --output=logs/05_syri_%j.log

# Boilerplate
set -euo pipefail
SIF=$(readlink -f images/sif/QC.sif)
TREM=$(readlink -f reference_data/lmultiflorum.tremona.primary.fa)
CIAO=$(readlink -f reference_data/ciao_unp.fasta)
OUTDIR="results/QC/syri"
T=${SLURM_CPUS_PER_TASK:-8}
BIND="/cluster/scratch"
CHR_REGEX="${CHR_REGEX:-^chr[0-9]+$}"
TEST_CHR="${TEST_CHR:-}"
run()    { singularity exec --bind "${BIND}" "${SIF}" "$@"; }
run_sh() { singularity exec --bind "${BIND}" "${SIF}" bash -c "$@"; }

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Extract chromosomes from both assemblies
if [[ ! -s tremona_chr.fa  ]]; then run seqkit grep -r -n -p "${CHR_REGEX}" "${TREM}" > tremona_chr.fa;  fi
if [[ ! -s brunharo_chr.fa ]]; then run seqkit grep -r -n -p "${CHR_REGEX}" "${CIAO}" > brunharo_chr.fa; fi

run samtools faidx tremona_chr.fa
run samtools faidx brunharo_chr.fa

if [[ ! -f aln.bam ]]; then
    if [[ ! -f aln.unsorted.bam ]]; then
        run_sh "minimap2 -ax asm5 -t $((T-2)) --eqx --secondary=no -L \
                  brunharo_chr.fa tremona_chr.fa \
                | samtools view -bS - > aln.unsorted.bam"
    fi
    run samtools sort -m 8G -@ ${T} -o aln.bam aln.unsorted.bam
    run samtools index -@ 4 aln.bam
    rm -f aln.unsorted.bam
fi

if [[ ! -f syri.out ]]; then
    run syri -c aln.bam -r brunharo_chr.fa -q tremona_chr.fa \
        -k -F B --nc 4
fi

# Plotsr genomes file
cat > genomes.txt <<'EOF'
#file	name	tags
brunharo_chr.fa	Brunharo	lw:1.5
tremona_chr.fa	Tremona	lw:1.5
EOF

# Whole-genome overview
run plotsr --sr syri.out --genomes genomes.txt -o syri_overview.png -d 300

# Per-chromosome plots
for chr in $(run seqkit fx2tab -n -i brunharo_chr.fa); do
    run plotsr --sr syri.out --genomes genomes.txt --chr "$chr" \
        -o "syri_${chr}.png" -d 300
done