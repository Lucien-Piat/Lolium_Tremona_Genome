import pandas as pd


def count_genes(gff):
    with open(gff) as fh:
        return sum(1 for l in fh if not l.startswith("#") and "\tgene\t" in l)


def measured_bp(fai):
    with open(fai) as fh:
        return sum(int(l.split("\t")[1]) for l in fh if l.strip())


def read_benchmark(path):
    df = pd.read_csv(path, sep="\t")
    return float(df.iloc[0]["s"]), float(df.iloc[0].get("max_rss", float("nan")))


rows = []
for panel, bench, gff, fai in (
        ("size", snakemake.input.bench, snakemake.input.gff, snakemake.input.fai),
        ("density", snakemake.input.dbench, snakemake.input.dgff,
         snakemake.input.dfai)):
    for b, g, f in zip(bench, gff, fai):
        secs, rss = read_benchmark(b)
        bp = measured_bp(f)
        rows.append({"panel": panel, "mb": bp / 1e6, "genes": count_genes(g),
                     "seconds": secs, "peak_rss_mb": rss})

df = pd.DataFrame(rows)
df["genes_per_mb"] = df.genes / df.mb
df["minutes"] = df.seconds / 60.0
df["peak_rss_gb"] = df.peak_rss_mb / 1024.0
df.to_csv(snakemake.output[0], sep="\t", index=False, float_format="%.4g")
