#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 01/06/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 3
#$ -t 1-22 -tc 5
#$ -N prune_extract
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load PLINK/2.00a2.3_x86_64

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/obesity_genotypes

mkdir tmp

# Get obesity-associated SNPs on this chromosome
awk -v pat="chr${SGE_TASK_ID}" '$0 ~ pat{print $4}' \
/well/lindgren/samvida/Resources/GWASCatalog/gwascat_obesity_associations_hg19.bed \
> tmp/chr${SGE_TASK_ID}_all_variants.txt

# Get genotype information at these SNPs 
plink2 \
--extract tmp/chr${SGE_TASK_ID}_all_variants.txt \
--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${SGE_TASK_ID}_v3.bgen ref-first \
--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--threads 3 \
--memory 15000 \
--make-pgen \
--out tmp/chr${SGE_TASK_ID}_all_variants

# Create LD-pruned list of independent variants (+/- 500kb, steps of 50, R2 < 0.5)
plink2 \
--pfile tmp/chr${SGE_TASK_ID}_all_variants \
--indep-pairwise 500 kb 0.5 \
--threads 3 \
--memory 15000 \
--out chr${SGE_TASK_ID}

# Make pgen files for obesity-associated variants on these chromosomes only
plink2 \
--extract chr${SGE_TASK_ID}.prune.in \
--pfile tmp/chr${SGE_TASK_ID}_all_variants \
--threads 3 \
--memory 15000 \
--make-pgen \
--out chr${SGE_TASK_ID}_pruned

# Remove genotype files for all variants (only retain pruned)
rm tmp/chr${SGE_TASK_ID}_all_variants.p*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
