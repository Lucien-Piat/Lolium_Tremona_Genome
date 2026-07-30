#!/bin/bash
set -euo pipefail

SIF=$(readlink -f images/sif/compleasm.sif)
ASM="./results/04_purgegrass/assembly.fa"
OUT="./results/04_purgegrass/compleasm_results/"
THREADS=8
run() { singularity exec "${SIF}" "$@"; }

mkdir -p "$OUT"

run compleasm run -a "$ASM" -o "$OUT" -l poales_odb12 -t "$THREADS"

