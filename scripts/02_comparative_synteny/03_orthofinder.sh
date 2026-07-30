#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
SIF="${ROOT}/images/sif/genome_analysis.sif"   # doit contenir orthofinder, sinon construis une image dediee
run(){ singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

PROT="${ROOT}/results/02_orthology/proteomes"
mkdir -p "${PROT}"
for n in tremona rabiosa sikem perenne paraquat brachypodium oryza; do
    cp "${ROOT}/results/02_synteny/${n}/proteins.faa" "${PROT}/${n}.faa"
done

run orthofinder -f "${PROT}" -t 4 -a 4
# resultat -> ${PROT}/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv