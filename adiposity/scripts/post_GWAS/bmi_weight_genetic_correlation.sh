#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 06/12/22

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N ldsc_bmi_weight
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/rg_h2

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

for SEX_STRATA in F M sex_comb; do
	for parameter in b1 k1 k1_k2 k1_k2_k3; do
		ldsc.py \
		--rg BMI_${SEX_STRATA}/tmp_${parameter}.sumstats.gz,Weight_${SEX_STRATA}/tmp_${parameter}.sumstats.gz \
		--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--out ${SEX_STRATA}_${parameter}_h2_rg_BMI_Weight
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
