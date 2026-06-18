#!/bin/bash

# Usage: bash scripts/te_annotation/09_gene_te.sh

set -euo pipefail
SIF=$(readlink -f images/sif/collapse_diag.sif)
BIND=$PWD
: "${TMPDIR:=/tmp}"
run() { singularity exec --bind "${BIND}" --bind "${TMPDIR}" "${SIF}" "$@"; }

OUTDIR=results/te_hite/gene_te
RMOUT=results/te_hite/tremona_TE.out.gz
GENEGFF=reference_data/lmultiflorum.tremona.gene_annotation.gff
FAI=reference_data/lmultiflorum.tremona.fa.fai
BM=results/dupclass/block_master.tsv
mkdir -p "${OUTDIR}"

cut -f1,2 "${FAI}" | sort -k1,1 > "${OUTDIR}/genome.txt"
read_out() { if [[ "${RMOUT}" == *.gz ]]; then zcat "${RMOUT}"; else cat "${RMOUT}"; fi; }

# TE BED: chrom start0 end family class div strand
read_out | awk 'NR>3 && NF>=11 {
        s=$6-1; e=$7; strand=($9=="C"?"-":"+")
        fam=$10; sub(/-(int|INT|LTR|I)$/, "", fam)
        printf "%s\t%d\t%d\t%s\t%s\t%s\t%s\n", $5, s, e, fam, $11, $2, strand
    }' | sort -k1,1 -k2,2n > "${OUTDIR}/te.bed"

# gene BED: chrom start0 end gene_id . strand
awk -F'\t' '$3=="gene" {
        match($9, /ID=[^;]+/); id=substr($9, RSTART+3, RLENGTH-3)
        printf "%s\t%d\t%d\t%s\t.\t%s\n", $1, $4-1, $5, id, $7
    }' "${GENEGFF}" | sort -k1,1 -k2,2n > "${OUTDIR}/genes.bed"

# nearest TE per gene
run bedtools closest -a "${OUTDIR}/genes.bed" -b "${OUTDIR}/te.bed" \
        -g "${OUTDIR}/genome.txt" -d -t first \
    | awk -F'\t' '{print $NF"\t"$11}' > "${OUTDIR}/gene_nearest_te.tsv"

# distance-binned class counts
awk -F'\t' '
    { d=$1; c=$2
      b=(d==0?"0_overlap": d<=1000?"1_0-1kb": d<=5000?"2_1-5kb": d<=20000?"3_5-20kb":"4_>20kb")
      n[b"\t"c]++ }
    END { for (k in n) print k"\t"n[k] }' \
    "${OUTDIR}/gene_nearest_te.tsv" | sort > "${OUTDIR}/gene_te_distance_bins.tsv"

# TE class composition inside an interval set vs genome background
awk '{bp[$5]+=$3-$2} END{for(c in bp) printf "genome\t%s\t%d\n", c, bp[c]}' \
    "${OUTDIR}/te.bed" > "${OUTDIR}/te_class_in_genome.tsv"

te_in_blocks() {
    local BLK=$1 LAB=$2
    [[ -s "${BLK}" ]] || { echo "skip ${LAB}: ${BLK} missing or empty"; return; }
    run bedtools intersect -a "${OUTDIR}/te.bed" -b "${BLK}" -wa \
        | awk -v lab="${LAB}" '{bp[$5]+=$3-$2} END{for(c in bp) printf "%s\t%s\t%d\n", lab, c, bp[c]}' \
        > "${OUTDIR}/te_class_in_${LAB}.tsv"
    echo "TE class bp inside ${LAB} -> ${OUTDIR}/te_class_in_${LAB}.tsv"
}

te_in_blocks results/popgen_mask/artifacts.bed artefacts

echo "Done."