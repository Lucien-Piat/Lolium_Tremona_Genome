#!/bin/bash
# Self synteny for all genomes
set -euo pipefail

ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"
LIB="${ROOT}/scripts/qc/circo_plot/lib"
OUTBASE="${ROOT}/results/synteny"
T=4

COLLAPSE_ISOFORMS=TRUE 

run() { singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

DATASETS=(
  "tremona|reference_data/lmultiflorum.tremona.fa|reference_data/lmultiflorum.tremona.gene_annotation.gff|tr"
  "rabiosa|reference_data/lmultiflorum.rabiosa.fa|results/annotation/tremona_to_lmultiflorum.rabiosa.gff.gz|ra"
  "paraquat|reference_data/lmultiflorum.paraquat.fasta|results/annotation/tremona_to_lmultiflorum.paraquat.gff.gz|pq"
  "perenne|reference_data/lmultiflorum.perenne.fa|results/annotation/tremona_to_lmultiflorum.perenne.gff.gz|pe"
  "brachypodium|reference_data/brachypodium.fna.gz|reference_data/brachypodium.gff.gz|bd"
  "oryza|reference_data/oryza.fna.gz|reference_data/oryza.gff.gz|os"
)

run_self_synteny() {
    local NAME="$1" GENOME_IN="$2" GFF_IN="$3" PREFIX="$4"
    local OUTDIR="${OUTBASE}/${NAME}"
    mkdir -p "${OUTDIR}"

    local GENOME
    if [[ "${GENOME_IN}" == *.gz ]]; then
        GENOME="${OUTDIR}/genome.fa"
        [ -s "${GENOME}" ] || zcat "${ROOT}/${GENOME_IN}" > "${GENOME}"
    else
        GENOME="${ROOT}/${GENOME_IN}"
    fi

    local GFF="${OUTDIR}/annotation.gff3"
    if [ ! -s "${GFF}" ]; then
        zcat -f "${ROOT}/${GFF_IN}" \
            | awk -F'\t' '
                BEGIN{OFS="\t"}
                /^#/ {print; next}
                # flag organelle sequences (NCBI region line carries genome=...)
                $3=="region" && $9 ~ /genome=(chloroplast|mitochondrion|plastid|apicoplast)/ { org[$1]=1; next }
                ($1 in org) {next}            # drop every feature on an organelle
                NF>=8 && $7=="?" {$7="+"}      # keep the trans-splice strand fix
                {print}
            ' > "${GFF}"
    fi

    cd "${OUTDIR}"

    # choose which gene models to keep
    if [ ! -s keep_ids.txt ]; then
        if [ "${COLLAPSE_ISOFORMS}" = true ]; then
            echo "[$(date)] keeping longest isoform per gene"
            run awk -F'\t' '
                function canon(s){ sub(/^rna-/,"",s); sub(/^gene-/,"",s); return s }
                $3=="mRNA" {
                    if ($9 ~ /valid_ORF=False/) next
                    id=""; par=""
                    if (match($9,/ID=([^;]+)/,a))     id=a[1]
                    if (match($9,/Parent=([^;]+)/,b)) par=b[1]
                    gene[id]=par; if (!(id in len)) len[id]=0
                }
                $3=="CDS" {
                    par=""
                    if (match($9,/Parent=([^;]+)/,b)) par=b[1]
                    len[par]+=($5-$4+1)
                }
                END {
                    for (m in gene){ g=gene[m]
                        if (!(g in best) || len[m]>blen[g]){ best[g]=m; blen[g]=len[m] } }
                    for (g in best) print canon(best[g])
                }' "${GFF}" > keep_ids.txt
        else
            echo "[$(date)] keeping all isoforms"
            run awk -F'\t' '
                function canon(s){ sub(/^rna-/,"",s); sub(/^gene-/,"",s); return s }
                $3=="mRNA" {
                    if ($9 ~ /valid_ORF=False/) next
                    if (match($9,/ID=([^;]+)/,a)) print canon(a[1])
                }' "${GFF}" > keep_ids.txt
        fi
    fi

    # proteins + MCScanX gff
    if [ ! -s proteins.raw.faa ]; then
        echo "[$(date)] gffread"
        run gffread -y proteins.raw.faa -g "${GENOME}" "${GFF}"
    fi
    if [ ! -s proteins.clean.faa ]; then
        echo "[$(date)] clean proteins"
        run python3 "${LIB}/clean_proteins.py" proteins.raw.faa proteins.clean.faa
    fi
    if [ ! -s proteins.faa ]; then
        echo "[$(date)] normalize headers + subset to kept models"
        run awk '
            function canon(s){ sub(/^rna-/,"",s); sub(/^gene-/,"",s); return s }
            NR==FNR { k[$0]=1; next }
            /^>/ { split(substr($0,2),p," "); id=canon(p[1]); keep=(id in k); if(keep) print ">" id; next }
            keep' keep_ids.txt proteins.clean.faa > proteins.faa
    fi
    if [ ! -s "${NAME}.gff" ]; then
        echo "[$(date)] MCScanX gff"
        run awk -F'\t' -v p="${PREFIX}" '
            function canon(s){ sub(/^rna-/,"",s); sub(/^gene-/,"",s); return s }
            NR==FNR { keep[$1]=1; next }
            $3=="mRNA" {
                if ($9 ~ /valid_ORF=False/) next
                id=""; if (match($9,/ID=([^;]+)/,a)) id=canon(a[1])
                if (!(id in keep)) next
                if (!($1 in seen)) { n++; seen[$1]=n }
                print p seen[$1] "\t" id "\t" $4 "\t" $5
            }' keep_ids.txt "${GFF}" > "${NAME}.gff"

    fi

    # DIAMOND self blast
    if [ ! -s proteins.db.dmnd ]; then
        echo "[$(date)] diamond makedb"
        run diamond makedb --in proteins.faa -d proteins.db --threads "${T}"
    fi
    if [ ! -s "${NAME}.blast" ]; then
        echo "[$(date)] diamond blastp"
        run diamond blastp \
            -q proteins.faa -d proteins.db -o "${NAME}.blast" \
            -e 1e-10 --outfmt 6 --max-target-seqs 5 --more-sensitive \
            --block-size 1.5 --index-chunks 2 --threads "${T}"
    fi

    # MCScanX + links
    if [ ! -s "${NAME}.collinearity" ]; then
        echo "[$(date)] MCScanX"
        run MCScanX -s 5 -e 1e-10 -m 25 "${OUTDIR}/${NAME}"
    fi
    if [ ! -s self_synteny_links.tsv ]; then
        echo "[$(date)] collinearity to links"
        run python3 "${LIB}/collinearity_to_links.py" \
            "${NAME}" "${PREFIX}" self_synteny_links.tsv
    fi

    echo "[$(date)] done ${NAME}"
    cd "${ROOT}"
}

mkdir -p "${OUTBASE}"
for d in "${DATASETS[@]}"; do
    IFS='|' read -r NAME GENOME GFF PREFIX <<< "${d}"
    run_self_synteny "${NAME}" "${GENOME}" "${GFF}" "${PREFIX}"
done

for d in "${DATASETS[@]}"; do
    IFS='|' read -r NAME _ _ _ <<< "${d}"
    COL="${OUTBASE}/${NAME}/${NAME}.collinearity"
    if [ -s "${COL}" ]; then
        N=$(grep -c "## Alignment" "${COL}" || true)
        echo "  ${NAME}: ${N} blocks"
    fi
done