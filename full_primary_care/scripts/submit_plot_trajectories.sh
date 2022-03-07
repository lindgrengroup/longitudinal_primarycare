#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 07/03/2022

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -N plot_gp_trajectories
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

# Input:
# 1. path to ids text file: (input format: text file with UKB eid and classification)
# 2. comma-separated list of biomarkers to plot ("all" to run all biomarkers in /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/code_lists/qcd_traits_available.txt)
# 3. name of log file
# 4. directory to where output plots should be stored

# If using bulk submission via R script:
# echo "passing covariates..."
# covars=${covars//|/,}
# echo ${covars}

logFile="/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/testing_plot_traj/sample_log.txt"
outPlotDir="/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/testing_plot_traj/plots/"

mkdir ${outPlotDir}

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/scripts/plot_trajectories.R \
--idFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/testing_plot_traj/sample_id_file_for_plot_trajectories.txt \
--biomarkers=all \
--logFile=${logFile} \
--outPlotDir=${outPlotDir}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
