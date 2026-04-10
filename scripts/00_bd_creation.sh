#!/bin/bash
set -euo pipefail

DBDIR="reference_data/taxdump"
mkdir -p "${DBDIR}"

echo "Downloading NCBI taxdump..."
wget -q --show-progress -O "${DBDIR}/new_taxdump.tar.gz" \
    https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar xzf "${DBDIR}/new_taxdump.tar.gz" -C "${DBDIR}"
rm -f "${DBDIR}/new_taxdump.tar.gz"

echo "Done. Upload to server:"
echo "  ${DBDIR}/"