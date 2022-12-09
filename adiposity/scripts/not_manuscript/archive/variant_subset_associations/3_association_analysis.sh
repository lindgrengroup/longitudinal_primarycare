#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 20/07/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 4
#$ -N plink_assocns_subset_geno
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

module load PLINK/2.00a2.3_x86_64

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/variant_subset_associations/
mkdir tmp_${STRATA}

# Run through each chromosome
for CHR in {1..22}; do 
	plink2 \
	--pfile /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/obesity_genotypes/chr${CHR}_pruned \
	--pheno /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/lmm_traits_${STRATA}.txt \
	--pheno-name lmm_intercepts,lmm_slopes_adj_int \
	--covar-name genotyping.array,UKB_assmt_centre,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
	--covar-variance-standardize \
	--threads 3 --memory 15000 \
	--glm hide-covar \
	--out tmp_${STRATA}/chr${CHR}
done

# Collate results from all chromosomes
# For lmm_intercepts
head -n 1 tmp_${STRATA}/chr1.lmm_intercepts.glm.linear > ${STRATA}_lmm_intercepts.txt
tail -n +2 -q tmp_${STRATA}/chr*.lmm_intercepts.glm.linear >> ${STRATA}_lmm_intercepts.txt

# For lmm_slopes
head -n 1 tmp_${STRATA}/chr1.lmm_slopes_adj_int.glm.linear > ${STRATA}_lmm_slopes_adj_int.txt
tail -n +2 -q tmp_${STRATA}/chr*.lmm_slopes_adj_int.glm.linear >> ${STRATA}_lmm_slopes_adj_int.txt

rm -r tmp_${STRATA}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
