#!/bin/bash

VERSION="0.2.7"
THREADS=4
OUTPUT_DIR="compleasm_results"
LINEAGE="poales"

mkdir -p "$OUTPUT_DIR"

for FASTA in *.fa.gz *.fasta.gz; do
    [ -e "$FASTA" ] || continue 
    
    BASENAME=$(basename "$FASTA" | sed -E 's/\.(fasta|fa)(\.gz)?$//')
    OUT_PATH="$OUTPUT_DIR/$BASENAME"
    SUMMARY_FILE="$OUT_PATH/summary.txt"
    
    echo "---------------------------------------------------"
    echo "Evaluating: $BASENAME"
    echo "---------------------------------------------------"
    
    if [ -f "$SUMMARY_FILE" ]; then
        echo "[SKIP] summary.txt found. $BASENAME is already processed."
        continue
    fi
    

    compleasm run \
        -a "$FASTA" \
        -o "$OUT_PATH" \
        -l "$LINEAGE" \
        -t "$THREADS"
        
    if [ $? -eq 0 ]; then
        echo "[CLEANUP] Removing copied ODB lineage from $OUT_PATH..."
        rm -rf "$OUT_PATH"/*_odb*
        echo "[SUCCESS] Finished $BASENAME."
    else
        echo "[ERROR] compleasm failed on $BASENAME. Check the logs."
    fi

done

echo "Done."