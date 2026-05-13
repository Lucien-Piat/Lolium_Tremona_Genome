#!/bin/bash
# Run compleasm locally
set -euo pipefail

ASM="./results/04c_purgegrass/assembly.fa"
OUT="./results/04c_purgegrass/compleasm_results/"
THREADS=8

mkdir -p "$OUT"

compleasm run -a $ASM -o "$OUT" -l poales -t "$THREADS"

