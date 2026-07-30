#!/bin/bash
# Ka/Ks (NG) inter/intra, tous les genomes sur la meme echelle
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/genome_analysis.sif"
LIB="${ROOT}/scripts/02_comparative_synteny"

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

DATASETS=(
  "tremona|reference_data/lmultiflorum.tremona.fa"
  "rabiosa|reference_data/lmultiflorum.rabiosa.fa"
  "sikem|reference_data/lmultiflorum.sikem.fasta"
  "paraquat|reference_data/lmultiflorum.paraquat.fasta"
  "perenne|reference_data/lmultiflorum.perenne.fa"
  "brachypodium|results/synteny/brachypodium/genome.fa"
  "oryza|results/synteny/oryza/genome.fa"
)

for d in "${DATASETS[@]}"; do
    IFS='|' read -r NAME GENOME <<< "${d}"
    SYN="${ROOT}/results/synteny/${NAME}"
    GFF="${SYN}/annotation.gff3"

    [ -s "${SYN}/cds.fa" ] && continue
    echo "[$(date)] CDS ${NAME}"

    # drop out-of-bounds coords 
    run awk -F'\t' '
        BEGIN{OFS="\t"}
        /^#/{print;next}
        $3=="region" && $9 ~ /genome=(chloroplast|mitochondrion|plastid|apicoplast)/ {org[$1]=1; next}
        ($1 in org){next}
        NF>=8 && $7=="?"{$7="+"}
        {print}
    ' "${GFF}" > "${SYN}/annotation.nuc.gff3"

    run gffread -x "${SYN}/cds.raw.fa" -g "${ROOT}/${GENOME}" "${SYN}/annotation.nuc.gff3"
    run awk '
        function canon(s){ sub(/^rna-/,"",s); sub(/^gene-/,"",s); return s }
        /^>/ { split(substr($0,2),p," "); print ">" canon(p[1]); next }
        { print }' "${SYN}/cds.raw.fa" > "${SYN}/cds.fa"
done
run python3 "${LIB}/05_synteny_ks_all.py"
echo "[$(date)] Done -> results/synteny/synteny_ks_all.pdf"