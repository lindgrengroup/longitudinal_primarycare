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

mkdir tmp

# Get obesity-associated SNPs on this chromosome
awk -v pat="chr${SGE_TASK_ID}" '$0 ~ pat{print $4}' \
/well/lindgren/samvida/Resources/GWASCatalog/gwascat_obesity_associations_hg19.bed \
> tmp/tmp_chr${SGE_TASK_ID}_variants.txt

# Get LD-pruned list of independent variants
plink2 \
--extract tmp/tmp_chr${SGE_TASK_ID}_variants.txt \
--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${SGE_TASK_ID}_v3.bgen ref-first \
--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--threads 3 \
--memory 15000 \
--make-pgen \
--out /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/obesity_genotypes/chr${SGE_TASK_ID}

# Get LD-pruned list of independent variants (+/- 500kb, steps of 50, R2 < 0.5)
plink2 \
--pfile /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/obesity_genotypes/chr${SGE_TASK_ID} \
--indep-pairwise 500 kb 0.5 \
--threads 3 \
--memory 15000 \
--out /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/obesity_genotypes/chr${SGE_TASK_ID}

rm tmp/tmp_chr${SGE_TASK_ID}_variants.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
