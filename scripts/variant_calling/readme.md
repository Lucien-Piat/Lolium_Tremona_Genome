

## Filters applied

**`cohort_allsites.vcf.gz`** No site filters. Contains every reference position, variant and invariant, unfiltered. The only sites missing are a small number
(~679 genome-wide) skipped by GATK for exceeding 50 alternate alleles.

**`cohort_snps_filtered.vcf.gz`** Derived from the all-sites file by:

```
bcftools view -v snps           # keep SNPs only 
bcftools filter -i 'QUAL>20 && INFO/QD>8'
```

- `QUAL > 20`  : phred-scaled site confidence above 20 (>99% confident a variant exists)
- `INFO/QD > 8`: quality normalised by depth above 8, removes variants whose high QUAL
  is only an artefact of very high coverage

Thresholds follow Stritt et al. 2022.
## Counts

- All-sites positions: 2,254,248,612
- Filtered SNPs: 113,268,915

