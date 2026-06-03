#!/bin/bash
# Ka/Ks (NG) inter/intra, tous les genomes sur la meme echelle
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
LIB="${ROOT}/scripts/qc/syntheny"

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

DATASETS=(
  "tremona|reference_data/lmultiflorum.tremona.fa"
  "rabiosa|reference_data/lmultiflorum.rabiosa.fa"
  "paraquat|reference_data/lmultiflorum.paraquat.fasta"
  "perenne|reference_data/lmultiflorum.perenne.fa"
  "brachypodium|results/synteny/brachypodium/genome.fa"
  "oryza|results/synteny/oryza/genome.fa"
)

for d in "${DATASETS[@]}"; do
    IFS='|' read -r NAME GENOME <<< "${d}"
    SYN="${ROOT}/results/synteny/${NAME}"
    GFF="${SYN}/annotation.gff3"

    if [ -s "${SYN}/cds.fa" ]; then
        continue
    fi
    echo "[$(date)] CDS ${NAME}"
    run gffread -x "${SYN}/cds.raw.fa" -g "${ROOT}/${GENOME}" "${GFF}"
    run awk '
        function canon(s){ sub(/^rna-/,"",s); sub(/^gene-/,"",s); return s }
        /^>/ { split(substr($0,2),p," "); print ">" canon(p[1]); next }
        { print }' "${SYN}/cds.raw.fa" > "${SYN}/cds.fa"
done

run python3 "${LIB}/synteny_ks_all.py"
echo "[$(date)] Done -> results/synteny/synteny_ks_all.pdf"