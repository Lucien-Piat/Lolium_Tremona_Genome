#!/bin/bash
set -euo pipefail

ASM="lmultiflorum.hap1.scaffolded.placed.fa"
OUT="compleasm"
THREADS=4

mkdir -p "$OUT"

compleasm run -a $ASM -o "$OUT" -l poales -t "$THREADS"

