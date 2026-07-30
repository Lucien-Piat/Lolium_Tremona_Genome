#!/bin/bash

INPUT_CSV="reads/paraquat.txt"

while IFS=',' read -r srr sample; do
    [ -z "$srr" ] && continue
    
    srr=$(echo "$srr" | tr -d '\r')
    sample=$(echo "$sample" | tr -d '\r')

    if [[ -f "${sample}_R1.fastq.gz" && -f "${sample}_R2.fastq.gz" ]]; then
        echo "Skipping $sample ($srr) - already downloaded and compressed."
        continue
    fi

    echo "Fetching $srr and renaming to $sample..."
    fasterq-dump --split-files --progress "$srr"

    if [[ -f "${srr}_1.fastq" && -f "${srr}_2.fastq" ]]; then
        mv "${srr}_1.fastq" "${sample}_R1.fastq"
        mv "${srr}_2.fastq" "${sample}_R2.fastq"
        echo "Success: Generated ${sample}_R1.fastq and ${sample}_R2.fastq"
        
        sbatch scripts/compress.sh "${sample}_R1.fastq"
        sbatch scripts/compress.sh "${sample}_R2.fastq"
        echo "Dispatched SLURM compression jobs for $sample."
        
    elif [[ -f "${srr}.fastq" ]]; then
        # Fallback for single-end datasets
        mv "${srr}.fastq" "${sample}.fastq"
        echo "Note: $srr appears to be single-end. Renamed to ${sample}.fastq"
        sbatch scripts/compress.sh "${sample}.fastq"
    else
        echo "Error: Failed to process $srr"
    fi
done < "$INPUT_CSV"

echo "Done."