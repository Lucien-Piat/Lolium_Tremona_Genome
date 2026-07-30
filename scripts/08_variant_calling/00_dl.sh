#!/bin/bash

INPUT_CSV="reads/paraquat.txt"

while IFS=',' read -r srr sample; do
    [ -z "$srr" ] && continue

    srr=$(echo "$srr" | tr -d '\r')
    sample=$(echo "$sample" | tr -d '\r')

    if [[ -f "${sample}_R1.fastq.gz" && -f "${sample}_R2.fastq.gz" ]]; then
        continue
    fi

    fasterq-dump --split-files --progress "$srr"

    if [[ -f "${srr}_1.fastq" && -f "${srr}_2.fastq" ]]; then
        mv "${srr}_1.fastq" "${sample}_R1.fastq"
        mv "${srr}_2.fastq" "${sample}_R2.fastq"
        sbatch scripts/08_variant_calling/compress.sh "${sample}_R1.fastq"
        sbatch scripts/08_variant_calling/compress.sh "${sample}_R2.fastq"
    else
        echo "Error: Failed to process $srr" >&2
    fi
done < "$INPUT_CSV"