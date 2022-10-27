#!/bin/bash
# Author: Samvida S. Venkatesh
# Date: 23/10/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N grep_snp_effects_all_strata
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/BOLT_results

awk 'NR>1{print $1}' \
/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/adipo_change/adipo_change_snps_replicate.txt \
> collated_lead_snps.txt

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb
do
	for CLUSTK in k1 k2 k3 k4
	do
		head -1 ${STRATA}_${CLUSTK}_final.txt > ${STRATA}_${CLUSTK}_all_leadsnp_assocns.txt
		grep -wFf collated_lead_snps.txt ${STRATA}_${CLUSTK}_final.txt \
		>> ${STRATA}_${CLUSTK}_all_leadsnp_assocns.txt
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
