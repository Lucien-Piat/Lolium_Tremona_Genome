import os
import subprocess
import sys

fasta = snakemake.input.fa
gff = snakemake.input.gff
out_bed = snakemake.output.bed
synteny = snakemake.output.synteny
outdir = os.path.dirname(out_bed)
os.makedirs(outdir, exist_ok=True)

SYNTENY_HEADER = ("block_id\tq_chr\tq_start\tq_end\tt_chr\tt_start\tt_end\t"
                  "median_ks\tage\n")
# Upstream failures that mean "no anchors", not "something broke". Both are
# reachable at the far end of the annotation-degradation series, where a
# sensitivity of zero is the measurement.
TOLERATED = ("'Sequence' is not in list",       # paralies/ks.py:142
             "Input file seems to be empty")    # diamond, no proteins


def record_empty():
    if not os.path.exists(out_bed):
        open(out_bed, "w").close()
    if not os.path.exists(synteny):
        with open(synteny, "w") as fh:
            fh.write(SYNTENY_HEADER)


def counts(path):
    genes = cds = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            genes += "\tgene\t" in line
            cds += "\tCDS\t" in line
    return genes, cds


n_genes, n_cds = counts(gff)
if n_genes == 0 or n_cds == 0:
    record_empty()
    sys.exit(0)

res = subprocess.run(
    ["/opt/conda/bin/paralies", "run", "--fasta", fasta, "--gff", gff,
     "--threshold", str(snakemake.params.threshold),
     "--threads", str(snakemake.threads), "--bed-only", "--outdir", outdir],
    capture_output=True, text=True)
out = res.stdout + res.stderr

if res.returncode == 0 or any(sig in out for sig in TOLERATED):
    record_empty()
    sys.exit(0)

sys.exit(f"ParaLies failed on {gff} (exit {res.returncode})\n"
         + "".join(out.splitlines(keepends=True)[-20:]))
