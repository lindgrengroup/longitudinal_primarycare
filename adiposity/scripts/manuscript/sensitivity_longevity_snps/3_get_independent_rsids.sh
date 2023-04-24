#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J get_longevity_rsids
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/logs/get_longevity_rsids-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load PLINK/2.00a2.3_x86_64

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity

# Get chr_pos for 1000G
awk -F '\t' 'BEGIN {OFS = FS} { NR>1 ; gsub(/chr/, ""); print $1"_"$3}' \
/well/lindgren/samvida/Resources/GWASCatalog/gwascat_longevity_associations_hg19.bed \
> tmp_longevity_chrpos.txt

# Only get independent SNPs (R2 < 0.1 with all other SNPs in +/- 500kb window)
# Based on LD patterns in 1000G EUR
echo EUR > tmp_EUR.txt
plink2 \
--bfile /well/lindgren/1000G_Phase3_chr_pos \
--keep-fam tmp_EUR.txt \
--extract tmp_longevity_chrpos.txt \
--indep-pairwise 1000kb 1 0.1 \
--out indep_longevity_chrpos \
--memory 15000 \
--threads 1 

# Get chromosome and position, sort by chr:pos
sort -k1,1 indep_longevity_chrpos.prune.in \
> tmp_indep_longevity_chrpos.txt

# Create chr_pos column in hg37 map file and sort by this
awk -F ' ' '{ NR>1; print $1"_"$2"\t"$3}' /well/lindgren/samvida/Resources/hg19/grch37.all.rsid.txt \
| sort -k1,1 > tmp_map.txt 

# Match chromosome and position to hg37 map file to get rsid
join -j 1 -o 1.1,2.2 tmp_indep_longevity_chrpos.txt tmp_map.txt \
> gwascat_indep_longevity_rsids.txt

rm tmp_*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
