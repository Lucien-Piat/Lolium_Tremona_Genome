#!/bin/bash
set -euo pipefail
GENOME=reference_data/lmultiflorum.tremona.fa
OUTDIR=results/te_hite/chunks
mkdir -p "${OUTDIR}"

cut -f1 "${GENOME}.fai" | while read -r name; do
    samtools faidx "${GENOME}" "${name}" > "${OUTDIR}/${name}.fa"
    echo "  ${name} -> ${OUTDIR}/${name}.fa"
done

echo "done."