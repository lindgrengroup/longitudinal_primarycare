#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 11/01/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1 
#$ -N power_boost
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/${STRATA_NAME}

# Extract SNPs with chisquare > 30 from GP GWAS file
awk -F '[[:space:]]+' 'NR>1 { if (($7/$8)^2 > 30) print $1 }' ${GP_FILENAME} \
> sig_snps_gp_gwas.txt

# Extract SNPs with chisquare > 30 from GIANT file
zcat ${META_FILENAME} | \
awk -F '[[:space:]]+' 'NR>1 { if (($7/$8)^2 > 30) print $3 }' - \
> sig_snps_giant_meta.txt

# Combine significant SNPs into single file
cat sig_snps_* > sig_snps.txt

# Extract sumstats for significant SNPs only
awk -F '\t' 'NR==FNR{a[$0]; next} FNR==1 || $1 in a' \
sig_snps.txt ${GP_FILENAME} > tmp_gp_gwas.txt

zcat ${META_FILENAME} | \
awk -F '[[:space:]]+' 'NR==FNR{a[$0]; next} FNR==1 || $3 in a' \
sig_snps.txt - > tmp_giant_meta.txt

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/8_calculate_power_boost.R \
--strata=${STRATA_NAME} \
--sample_size_gp=${GP_N} 

rm tmp_*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
