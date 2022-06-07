#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 01/06/21

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 1
#$ -t 1-12 -tc 3
#$ -N metabo_assocns_lmms

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load PLINK/2.00a2.3_x86_64

##########################################################################################

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/strata_filenames.txt`
mkdir /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/PLINK_results/${STRATA_NAME}

for c in {1..22}
do
	plink2 \
	--extract /well/lindgren/UKBIOBANK/samvida/full_primary_care/GWAS/metabo_endo_snps_passed_QC_211102/passed_QC.txt \
	--pfile /well/lindgren/UKBIOBANK/samvida/full_primary_care/GWAS/metabo_endo_snps_passed_QC_211102/chr${c} \
	--chr ${c} \
	--covar /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/sample_qc/${STRATA_NAME}_ids_passed_qc.txt \
	--glm hide-covar \
	--covar-variance-standardize \
	--threads 1 \
	--memory 15000 \
	--pheno /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/all_models_${STRATA_NAME}.txt \
	--out /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/PLINK_results/${STRATA_NAME}/assoc_chr${c}
done

cd /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/PLINK_results/${STRATA_NAME}

TRAITS=('lmm_intercepts' 'lmm_slopes_adj_baseline')
 
for trait in "${TRAITS[@]}"
do
	awk '(NR == 1) || (FNR > 1)' assoc_chr*.${trait}.glm.linear > ../${STRATA_NAME${trait}_assoc_results.txt
	rm assoc_chr*.${trait}.glm.linear
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0