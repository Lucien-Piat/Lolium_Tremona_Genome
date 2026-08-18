import pandas as pd

frames = [pd.read_csv(p, sep="\t") for p in snakemake.input]
df = pd.concat(frames, ignore_index=True)
df.to_csv(snakemake.output[0], sep="\t", index=False, float_format="%.6g")
