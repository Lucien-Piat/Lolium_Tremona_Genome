#!/usr/bin/env bash
#SBATCH --job-name=popgen_mask
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --time=03:00:00
#SBATCH --output=logs/mask_%j.log

module load samtools
module load bedtools
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snakemake

set -euo pipefail

REF=mask/tremona_TE.softmasked.fa
TE_GFF=mask/tremona_TE.gff3
ART_BED=mask/masked_intervals.bed
OUTDIR=results/popgen_mask
READLEN=150
THREADS=${SLURM_CPUS_PER_TASK:-4}

mkdir -p "$OUTDIR"
[ -f "${REF}.fai" ] || samtools faidx "$REF"
cut -f1,2 "${REF}.fai" | sort -k1,1 > "$OUTDIR/genome.txt"

awk -F'\t' 'BEGIN{OFS="\t"} /^##FASTA/{exit} /^#/{next}
     NF>=8 && $4 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/ {print $1,$4-1,$5}' "$TE_GFF" \
  | sort -k1,1 -k2,2n > "$OUTDIR/te.bed"

# artefacts
cut -f1-3 "$ART_BED" | sort -k1,1 -k2,2n > "$OUTDIR/artifacts.bed"

# mappability
[ -d "$OUTDIR/genmap_idx" ] || genmap index -F "$REF" -I "$OUTDIR/genmap_idx"
genmap map -K "$READLEN" -E 0 -T "$THREADS" -I "$OUTDIR/genmap_idx" -O "$OUTDIR/genmap" -bg
awk -F'\t' '$4 < 1 {print $1"\t"$2"\t"$3}' "$OUTDIR/genmap.bedgraph" \
  | sort -k1,1 -k2,2n | bedtools merge -i - > "$OUTDIR/lowmap.bed"

cat "$OUTDIR"/{te,artifacts,lowmap}.bed | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$OUTDIR/exclude.bed"
bedtools sort -i "$OUTDIR/exclude.bed" -g "$OUTDIR/genome.txt" \
  | bedtools complement -i - -g "$OUTDIR/genome.txt" > "$OUTDIR/accessible.bed"

G=$(awk '{s+=$2} END{print s}' "$OUTDIR/genome.txt")
{
  printf "layer\tn_regions\tbp\tpct_genome\n"
  for f in te artifacts lowmap exclude accessible; do
    bedtools merge -i "$OUTDIR/$f.bed" \
      | awk -v lyr="$f" -v g="$G" '{n++; bp+=$3-$2} END{printf "%s\t%d\t%d\t%.3f\n", lyr, n+0, bp+0, 100*(bp+0)/g}'
  done
} > "$OUTDIR/mask_recap.tsv"

column -t "$OUTDIR/mask_recap.tsv"