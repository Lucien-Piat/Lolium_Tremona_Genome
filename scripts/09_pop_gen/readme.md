## 1. Masque d'accessibilite (results/mask/)

- **Mappabilite** : retirer les multimap (k=150, E=2, m >= 0.5)
  - retire : ~1 281 Mb (62.2 %)
  - retenu : 777 Mb (37.8 %)
- **TE** : retirer les TE
  - retire : 1 358 Mb (66.0 %)
  - retenu : 700 Mb (34.0 %)
- **Accessible** : intersection mappable ET non-TE
  - retire : 1 657 Mb (80.5 %)
  - retenu : **401 Mb (19.5 %)**

## 2. VCF filtres (results/filtered_vcf/)

### allsites.hc.vcf.gz -- pixy (Pi, Tajima's D)

- **Masque** : regions accessibles -> `-R accessible.bed`
- **Genotype DP** : appels faible profondeur -> `DP <= 10 -> ./.`
- **Genotype GQ** : appels faible qualite -> `GQ <= 30 -> ./.`

### snps.ld.vcf.gz -- PCA, structure sNMF, FIS

- **Depart** : SNP du joint calling (`cohort_snps_filtered`)
- **Masque + biallel** : `-R mask, -m2 -M2 -v snps`
- **MAF** : retirer variants rares -> `MAF >= 0.05` -> 18 050 507 SNP
- **Pruning LD** : retirer redondance -> `plink2 r2 > 0.1, 50 kb` -> 134 262 SNP
- **(option) Thin** : tirage aleatoire -> `plink2 --thin-count 50000` -> ~50 000 SNP

Evolution : 18 050 507 -> (r2>0.4 : 1 509 718) -> (r2>0.1 : 134 262) 
