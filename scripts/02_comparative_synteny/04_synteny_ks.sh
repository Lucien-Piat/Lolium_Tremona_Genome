#!/bin/bash
# Ka/Ks (NG) + overlay BUSCO, Tremona uniquement
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
LIB="${ROOT}/scripts/02_comparative_synteny"
SYN="${ROOT}/results/synteny/tremona_purged"
GENOME="${ROOT}/reference_data/lmultiflorum.tremona.fa"
GFF="${SYN}/annotation.gff3"
BUSCO="${ROOT}/reference_data/lmultiflorum.tremona_full_table_busco_format.tsv"

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

if [ ! -s "${SYN}/cds.fa" ]; then
    echo "[$(date)] Extracting CDS"
    run gffread -x "${SYN}/cds.raw.fa" -g "${GENOME}" "${GFF}"
    run awk '
        function canon(s){ sub(/^rna-/,"",s); sub(/^gene-/,"",s); return s }
        /^>/ { split(substr($0,2),p," "); print ">" canon(p[1]); next }
        { print }' "${SYN}/cds.raw.fa" > "${SYN}/cds.fa"
fi

N=$(comm -12 \
    <(awk -F'\t' '$3=="mRNA"{print $1}' "${GFF}" | sort -u) \
    <(awk -F'\t' '!/^#/{print $3}' "${BUSCO}" | sort -u) | wc -l)
echo "[$(date)] chromosomes communs GFF/BUSCO : ${N}"
[ "${N}" -eq 0 ] && echo "ATTENTION: overlay BUSCO vide, GFF et BUSCO sur assemblages differents." >&2

run python3 "${LIB}/04_synteny_ks.py" --keep-tmp
echo "[$(date)] Done -> ${SYN}/synteny_ks.pdf"