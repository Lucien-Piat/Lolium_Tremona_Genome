# Lolium multiflorum Tremona Genome

Chromosome-scale genome assembly, duplication curation, annotation, and
population genomics of the Swiss *L. multiflorum* ecotype Tremona. 

User friendly duplication curation pipeline is aviable here : [Paralies](https://github.com/Lucien-Piat/ParaLies)
## Running

Every script runs a tool through a Singularity/Apptainer container:
`images/def/*.def` are the recipes, built into `images/sif/*.sif`. Run scripts
from the repo root, in numeric order within each `scripts/NN_*/` folder
(`sbatch` for SLURM jobs, `bash` otherwise). `scripts/draft_*` and
`draft_*` files were exploratory and are not part of the paper.

## Steps

1. **assembly** — HiFi reads -> hifiasm -> purge haplotigs -> organelles (Oatk) -> RagTag scaffold
2. **comparative_synteny** — self-synteny + Ks across the 7-genome panel (Fig. 1)
3. **duplication_classification** — Ks/depth/shared-polymorphism artefact calling (Fig. 2a-d)
4. **surgery** — excise artefacts, fix orientation, before/after check (Fig. 2e-g, Fig. S5)
5. **assembly_qc** — Merqury QV (Fig. S3) + Circos genome landscape (Fig. 3)
6. **gene_annotation** — Liftoff gene models from GULF + Kyuss
7. **te_annotation** — HiTE de novo + RepeatMasker, repeat landscape (Fig. S4)
8. **variant_calling** — fastp -> minimap2 -> GATK, per-sample to joint genotyping
9. **pop_gen** — accessibility mask, PCA/sNMF/NJ tree, FIS, diversity, IBD (Fig. 4)
