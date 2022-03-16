#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 07/03/2022

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -N plot_biomarker_disease_trajectories
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

# Input:
# 1. comma-separated list of diagnoses to plot ("all" to run all diseases in /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt)
# 2. name of biomarker to plot
# 3. name of log file
# 4. directory where output plots should be stored

# If using bulk submission via R script:
echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

mkdir ${outPlotDir}

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/scripts/plot_trajectories_with_diagnosis_info.R \
--diagnoses=${diagnoses} \
--biomarker=${biomarker} \
--logFile=${logFile} \
--outPlotDir=${outPlotDir}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
