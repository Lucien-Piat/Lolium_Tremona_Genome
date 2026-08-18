# ParaLies simulation benchmark

Snakemake workflow behind the two simulation figures: duplications of known
divergence, size and age are injected into a pre-purged *Brachypodium
distachyon* assembly and recovered with ParaLies.
```
snakemake --cores 4
```
Outputs `results/figure_summary.pdf` and `results/figure_bp_errors.pdf`, plus
the tables behind them (`results/blocks_all.tsv`, `results/bp_all.tsv`,
`results/scaling.tsv`).

`config.yaml` states the design: the injection grids, the heterozygosity
levels, the degradation series and the scaling ladders. Every condition is run
`replicates` times over independent injections.
