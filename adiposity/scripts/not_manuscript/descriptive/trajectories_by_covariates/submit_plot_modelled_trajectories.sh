#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 07/03/2022

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -N plot_modelled_trajectories
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

# Input:
# 1. path to ids text file: (input format: text file with UKB eid and classification)
# 2. comma-separated list of biomarkers to plot ("all" to run BMI, weight, waist circumference, and WHR)
# 3. name of log file
# 4. directory to where output plots should be stored

# If using bulk submission via R script:
echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/plot_modelled_trajectories.R \
--idFile=${idFile} \
--biomarkers=${phenotype} \
--outPrefix=${outPrefix} 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
