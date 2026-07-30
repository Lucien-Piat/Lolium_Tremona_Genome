#!/bin/bash
# Run compleasm locally
set -euo pipefail

ASM="./results/04_purgegrass/assembly.fa"
OUT="./results/04_purgegrass/compleasm_results/"
THREADS=8

mkdir -p "$OUT"

compleasm run -a $ASM -o "$OUT" -l poales -t "$THREADS"

