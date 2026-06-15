#!/bin/bash

# Usage : bash scripts/te_annotation/07_class_table.sh

set -euo pipefail
OUTDIR=results/te_hite/annotation/merged
OUT="${OUTDIR}/tremona_TE.out"
GENOME=reference_data/lmultiflorum.tremona.primary.fa
TABLE="${OUTDIR}/tremona_TE.class_table.tsv"
TMP="${OUTDIR}/.cls.tmp"

GSIZE=$(awk '{s+=$2} END{print s}' "${GENOME}.fai")

# Somme bp et nombre de fragments par classe/famille colonne 11 
awk -v g="${GSIZE}" '
    NR>3 && NF>=11 { len=$7-$6+1; bp[$11]+=len; n[$11]++; tot+=len; ntot++ }
    END {
        for (c in bp) printf "C\t%s\t%d\t%d\t%.2f\n", c, n[c], bp[c], 100*bp[c]/g
        printf "T\t-\t%d\t%d\t%.2f\n", ntot, tot, 100*tot/g
    }' "${OUT}" > "${TMP}"

# Table triee par bp decroissant
{
    printf "class/family\tfragments\tbp\t%%genome\n"
    grep '^C' "${TMP}" | cut -f2- | sort -t$'\t' -k3,3 -rn
} > "${TABLE}"

column -t -s$'\t' "${TABLE}"
echo
grep '^T' "${TMP}" | awk -F'\t' \
    '{printf "TOTAL interspersed : %d fragments, %d bp, %.2f%% du genome\n", $3, $4, $5}'
rm -f "${TMP}"
echo "Table -> ${TABLE}"