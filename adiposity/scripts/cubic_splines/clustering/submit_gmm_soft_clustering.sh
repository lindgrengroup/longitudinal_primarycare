#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -N mclust_plots
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/gmm_soft_clustering.R \
--strata=$STRATA \
--nClust=$NCLUST 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
