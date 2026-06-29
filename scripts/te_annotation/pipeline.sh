#!/bin/bash
# scripts/te_annotation/run_te_pipeline.sh
# Usage: bash scripts/te_annotation/run_te_pipeline.sh

set -euo pipefail

BASEDIR=results/te_hite
RMOUT="${BASEDIR}/tremona_TE.out"
GENOME=reference_data/lmultiflorum.tremona.fa
FAI="${GENOME}.fai"
GENEGFF=reference_data/lmultiflorum.tremona.gene_annotation.gff
GENE_TE_DIR="${BASEDIR}/gene_te"
FIG_DIR="${BASEDIR}/figures"
mkdir -p "${GENE_TE_DIR}" "${FIG_DIR}"
BIND=$PWD
: "${TMPDIR:=/tmp}"
SIF=$(readlink -f images/sif/collapse_diag.sif)
PLOT_SIF=$(readlink -f images/sif/plotting.sif 2>/dev/null || echo "")
run()  { singularity exec --bind "${BIND}" --bind "${TMPDIR}" "${SIF}" "$@"; }
plot() {
    if [[ -n "${PLOT_SIF}" && -f "${PLOT_SIF}" ]]; then
        singularity exec --bind "${BIND}" --bind "${TMPDIR}" "${PLOT_SIF}" "$@"
    else
        "$@"
    fi
}

read_out() { if [[ "${RMOUT}" == *.gz ]]; then zcat "${RMOUT}"; else cat "${RMOUT}"; fi; }
GSIZE=$(awk '{s+=$2} END{print s}' "${FAI}")

# class composition table (class/family, fragments, bp, % genome)
CLASSTAB="${BASEDIR}/tremona_TE.class_table.tsv"
TMP="${BASEDIR}/.cls.tmp"

read_out | awk -v g="${GSIZE}" '
    NR>3 && NF>=11 { len=$7-$6+1; bp[$11]+=len; n[$11]++; tot+=len; ntot++ }
    END {
        for (c in bp) printf "C\t%s\t%d\t%d\t%.4f\n", c, n[c], bp[c], 100*bp[c]/g
        printf "T\t-\t%d\t%d\t%.4f\n", ntot, tot, 100*tot/g
    }' > "${TMP}"

{
    printf "klass\tfragments\tbp\tpct_genome\n"
    awk -F'\t' '$1=="C"' "${TMP}" | cut -f2- | sort -t$'\t' -k3,3rn
} > "${CLASSTAB}"

echo
echo "Repeat composition by class:"
{
    printf "class/family\tfragments\tlength_Mb\t%%genome\n"
    tail -n +2 "${CLASSTAB}" \
        | awk -F'\t' '{printf "%s\t%d\t%.2f\t%.2f\n", $1, $2, $3/1e6, $4}'
} | column -t -s$'\t'
awk -F'\t' '$1=="T" {printf "\nTOTAL interspersed: %d fragments, %.2f Mb, %.2f%% of genome\n", $3, $4/1e6, $5}' "${TMP}"
rm -f "${TMP}"
echo "Class table -> ${CLASSTAB}"

# TE bed, gene bed, chromosome sizes
cut -f1,2 "${FAI}" | sort -k1,1 > "${GENE_TE_DIR}/genome.txt"

read_out | awk 'NR>3 && NF>=11 {
        s=$6-1; e=$7; strand=($9=="C"?"-":"+")
        fam=$10; sub(/-(int|INT|LTR|I)$/, "", fam)
        printf "%s\t%d\t%d\t%s\t%s\t%s\t%s\n", $5, s, e, fam, $11, $2, strand
    }' | sort -k1,1 -k2,2n > "${GENE_TE_DIR}/te.bed"

awk -F'\t' '$3=="gene" {
        match($9, /ID=[^;]+/); id=substr($9, RSTART+3, RLENGTH-3)
        printf "%s\t%d\t%d\t%s\t.\t%s\n", $1, $4-1, $5, id, $7
    }' "${GENEGFF}" | sort -k1,1 -k2,2n > "${GENE_TE_DIR}/genes.bed"

# nearest gene distance
DIST="${GENE_TE_DIR}/te_gene_distance.tsv"
run bedtools closest -a "${GENE_TE_DIR}/te.bed" -b "${GENE_TE_DIR}/genes.bed" \
        -g "${GENE_TE_DIR}/genome.txt" -d -t first \
    | awk -F'\t' 'BEGIN{OFS="\t"} {len=$3-$2; print $NF, $5, len}' \
    > "${DIST}"
echo "TE->gene distance -> ${DIST}"

PART="${BASEDIR}/genome_partition.tsv"
CODING_BP=$(awk -F'\t' '$3=="CDS" {printf "%s\t%d\t%d\n", $1, $4-1, $5}' "${GENEGFF}" \
    | sort -k1,1 -k2,2n \
    | run bedtools merge -i - \
    | run bedtools subtract -a - -b "${GENE_TE_DIR}/te.bed" \
    | awk '{s+=$3-$2} END{print s+0}')

{
    printf "component\tbp\n"
    printf "genome_total\t%d\n" "${GSIZE}"
    printf "coding_nonTE\t%d\n" "${CODING_BP}"
} > "${PART}"
echo "Genome partition -> ${PART}  (genome=${GSIZE} bp, coding_nonTE=${CODING_BP} bp)"


plot python3 scripts/te_annotation/plot_te_composition.py
plot python3 scripts/te_annotation/plot_chr1_te_families.py
plot python3 scripts/te_annotation/plot_gene_te.py

echo "Done."