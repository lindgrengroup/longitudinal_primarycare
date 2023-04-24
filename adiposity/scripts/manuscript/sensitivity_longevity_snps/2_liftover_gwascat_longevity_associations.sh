#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J liftover_gwascat_longevity_associations
#SBATCH -o /well/lindgren/samvida/Resources/GWASCatalog/liftover_gwascat_longevity_associations-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
module load liftOver/20210519

cd /well/lindgren/samvida/Resources/GWASCatalog

# Extract variants associated with hormones
head -n 1 all_associations_230327.txt > gwascat_longevity_associations.txt
grep "longevity" all_associations_230327.txt >> gwascat_longevity_associations.txt

# Create bed file to liftover
Rscript scripts/1_make_gwascat_bed.R

liftOver gwascat_longevity_associations_hg38.bed \
/well/lindgren/samvida/Resources/hg38ToHg19.over.chain.gz \
gwascat_longevity_associations_hg19.bed \
gwascat_longevity_associations_unlifted.bed

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
