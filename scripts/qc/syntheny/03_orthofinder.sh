#!/bin/bash
set -euo pipefail
ROOT=$(pwd)
SIF="${ROOT}/images/sif/QC.sif"   # doit contenir orthofinder, sinon construis une image dediee
run(){ singularity exec --bind "${ROOT}":"${ROOT}" "${SIF}" "$@"; }

PROT="${ROOT}/results/orthology/proteomes"
mkdir -p "${PROT}"
for n in tremona rabiosa perenne paraquat brachypodium; do
    cp "${ROOT}/results/synteny/${n}/proteins.faa" "${PROT}/${n}.faa"
done

run orthofinder -f "${PROT}" -t 4 -a 4
# resultat -> ${PROT}/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv