#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 20/05/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -t 1-18 -tc 6
#$ -N Baseline_flexmix
#$ -j y

echo "########################################################"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started date: "`date`
echo "##########################################################"

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" strata_filenames.txt`

module load R-bundle-Bioconductor/3.9-foss-2019a-R-3.6.0

Rscript 5_baseline_models.R BMI $STRATA_NAME

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0