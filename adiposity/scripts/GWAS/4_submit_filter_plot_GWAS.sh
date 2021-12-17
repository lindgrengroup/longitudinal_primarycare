#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -hold_jid lmm_intercepts_GWAS
#$ -t 1-12 -tc 4 
#$ -N adipo_intercept_GWAS_plot
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" strata_filenames.txt`

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
Rscript 4_filter_GWAS_results.R $STRATA_NAME

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
