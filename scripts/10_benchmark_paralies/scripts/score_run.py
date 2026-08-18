import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

THRESHOLD = float(snakemake.params.threshold)
# Within this of the threshold the estimator's own error decides the call,
# so those blocks are counted apart from genuine over-purging.
BOUNDARY = 0.025

truth = common.read_truth(snakemake.input.truth)
bed = common.read_bed(snakemake.input.bed)
syn = getattr(snakemake.input, "synteny", None)
estimated = common.match_synteny(truth, syn) if syn else [float("nan")] * len(truth)

rep = str(snakemake.params.rep)
condition = str(snakemake.params.condition)
experiment = str(snakemake.params.experiment)


rows = []
for b, est in zip(truth.itertuples(), estimated):
    rows.append({
        "experiment": experiment, "rep": rep, "condition": condition,
        "injected_ks": b.injected_ks, "type": b.type, "n_genes": b.n_genes,
        "detected": int(common.block_detected(bed, b)),
        "estimated_ks": est,


        "copy_bp": min(b.src_end - b.src_start, b.dup_end - b.dup_start),
    })
blocks = common.pd.DataFrame(rows)
blocks.to_csv(snakemake.output.blocks, sep="\t", index=False,
              float_format="%.5g")


removed = sum(e - s for ivs in bed.values() for s, e in ivs)
correct = boundary = wrong = should = 0
for b in truth.itertuples():
    cov = (common.covered(bed, b.src_chrom, b.src_start, b.src_end)
           + common.covered(bed, b.dup_chrom, b.dup_start, b.dup_end))
    # Ground truth is the injected Ks, not the block label: the ks ladder
    # types everything "artefact" but spans 0 to 2.0.
    if b.type == "paralogue":
        wrong += cov
    elif b.injected_ks < THRESHOLD:
        correct += cov

        should += min(b.src_end - b.src_start, b.dup_end - b.dup_start)
    elif b.injected_ks <= THRESHOLD + BOUNDARY:
        boundary += cov
    else:
        wrong += cov

asm = sum(int(l.split("\t")[1]) for l in open(snakemake.input.fai) if l.strip())
injected = sum((b.src_end - b.src_start) + (b.dup_end - b.dup_start)
               for b in truth.itertuples())

common.pd.DataFrame([{
    "experiment": experiment, "rep": rep, "condition": condition,
    "removed_bp": removed, "correct_bp": correct, "boundary_bp": boundary,
    "wrong_bp": wrong, "background_bp": removed - correct - boundary - wrong,
    "should_remove_bp": should, "keep_bp": asm - injected, "assembly_bp": asm,
}]).to_csv(snakemake.output.bp, sep="\t", index=False, float_format="%.6g")

