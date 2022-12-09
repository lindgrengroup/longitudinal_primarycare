#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 07/02/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N polmm_step1
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_logs/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_logs/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/POLMM_helpers/1_POLMM_step1_wrapper.R \
	--strata=${STRATA} \
	--covars=baseline_age,age_sq,UKB_assmt_centre,genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
	--sampleIDColinphenoFile=eid \
	--ordinalColinphenoFile=clust \
	--outputdir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_results/ 
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
