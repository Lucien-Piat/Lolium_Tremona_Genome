Prepare the reads 
sbatch scripts/00_bam2fastq

Verify ploidy and read error
sbatch scripts/02_meryl_genomescope.sh

Run the assembly
sbatch scripts/03_hifiasm.sh

from https://github.com/Yutang-ETH/PhaseGrass/wiki/Step-2:-generating-a-chromosome%E2%80%90level-unphased-haploid-assembly

Purge the duplicated tigs
sbatch scripts/04c_pg_01_prep.sh

All to all alignement

hap ~27, dip ~52
low = 5 
mid = 40
high = 100

sbatch scripts/04c_pg_02_purgehap.sh
Busco with poales 
sbatch scripts/04c_pg_03_busco.sh



bash scripts/05_mitohifi.sh
