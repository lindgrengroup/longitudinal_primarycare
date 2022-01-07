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

cd /well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS

cat finemap_filenames.txt | while read STRATA_NAME || [[ -n $STRATA_NAME ]];
do
	# Get sample IDs to retain for finemapping
	awk 'NR > 1 {print $1}' \
	./sample_qc/${STRATA_NAME}_ids_passed_qc.txt \
	> ./sample_qc/${STRATA_NAME}_ids.incl
	
	# Print number of samples 
	wc -l ./sample_qc/${STRATA_NAME}_ids.incl
done

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/pipeline

source activate pipeline

cat /well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/finemap_filenames.txt | while read STRATA_NAME || [[ -n $STRATA_NAME ]];
do
	pipeline.r --outname=${STRATA_NAME} \
	--outdir=/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/post_GWAS/${STRATA_NAME} \
	--input-json=/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/post_GWAS/input_jsons/${STRATA_NAME}_inputs.json \
	--job-dependencies-json=./json/job_dependencies.json \
	--sumstats=/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/BOLT_filtered/${STRATA_NAME}.txt.gz
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
