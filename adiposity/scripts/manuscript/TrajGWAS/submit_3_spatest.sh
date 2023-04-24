#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 01/06/2022

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -t 1-22 
#$ -N hormone_trajgwas
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module purge
module load Julia/1.6.1-linux-x86_64

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

# Strata called by submission script
# chromosome from task ID
CHR=${SGE_TASK_ID}

julia \
/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/TrajGWAS/scripts/3_spatest.jl \
${STRATA} ${CHR}  

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
