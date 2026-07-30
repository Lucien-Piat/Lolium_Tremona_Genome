#!/bin/bash

ASSEMBLIES=$(cut -f1 datasets.tsv | tr '\n' ' ')

for file in $ASSEMBLIES; do
    if [ ! -f "$file" ]; then
        echo "Error: Assembly file not found -> $file" >&2
        exit 1
    fi
done

LABELS=$(cut -f2 datasets.tsv | paste -sd, -)

quast \
    $ASSEMBLIES \
    -o results/quast_output \
    --threads 4 \
    --large \
    --split-scaffolds \
    --no-snps \
    --labels "$LABELS"