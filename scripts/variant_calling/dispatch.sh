#!/bin/bash

# Usage: ./dispatch_pipeline.sh <SAMPLE_ID> [MAIN_JOB_ID]

SAMPLE=$1
MAIN_JOB_ID=$2

echo "Submitting pipeline for sample: $SAMPLE"

JOB1_ARGS=("--parsable")

if [ -n "$MAIN_JOB_ID" ]; then
    JOB1_ARGS+=("--dependency=afterok:$MAIN_JOB_ID")
fi

JOB1=$(sbatch "${JOB1_ARGS[@]}" ./scripts/01_align.sh "$SAMPLE")
echo "  -> Alignment submitted (Job ID: $JOB1)"

JOB2=$(sbatch --parsable --dependency=afterok:"$JOB1" ./scripts/02_markdup.sh "$SAMPLE")
echo "  -> MarkDup submitted (Job ID: $JOB2)"

JOB3=$(sbatch --parsable --dependency=afterok:"$JOB2" ./scripts/03_haplotype_caller.sh "$SAMPLE")
echo "  -> HaplotypeCaller submitted (Job ID: $JOB3)"