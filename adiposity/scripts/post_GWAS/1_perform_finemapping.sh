#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N submit_finemapping
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/pipeline

source activate pipeline

cat /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/strata_filenames.txt | while read STRATA || [[ -n $STRATA ]]; 
do
	# Quantitative traits from BOLT output
	for PARAMETER in "lmm_intercepts" "lmm_slopes_adj_baseline" "cspline_intercepts"
	do
		pipeline.r --outname=${STRATA}_${PARAMETER} \
		--outdir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/post_GWAS/${STRATA}/${PARAMETER} \
		--input-json=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/post_GWAS/0_input_jsons/${STRATA}_inputs.json \
		--job-dependencies-json=/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/pipeline/json/job_dependencies.json \
		--sumstats=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_results/${STRATA}/${STRATA}_${PARAMETER}_final.txt.gz
	done
	# Cluster traits from SAIGE output
	for KI in {1..6}
	do
		pipeline.r --outname=${STRATA}_k${KI} \
		--outdir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/post_GWAS/${STRATA}/k${KI} \
		--input-json=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/post_GWAS/0_input_jsons/${STRATA}_inputs.json \
		--job-dependencies-json=/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/pipeline/json/job_dependencies.json \
		--sumstats=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/SAIGE_results/${STRATA}/${STRATA}_k${KI}_final.txt.gz
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
