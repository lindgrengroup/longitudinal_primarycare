#!/bin/bash

# Author: George Nicholson
# Date: 16/02/23

# Specify working directory when executing script
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/georgenicholson/github_repos/longitudinal_primarycare/biobank_completion
# Specify which project your job belongs to
#$ -P lindgren.prjc 
# Specify which queue to run a job in. 
#$ -q short.qc
# Number of slots
#$ -pe shmem 1
# How many tasks in array (first id)-(last id):(step size)
#$ -t 1-100:1
# Number of tasks eligible for concurrent execution
#$ -tc 100
# Name of the array job
#$ -N bash_03_optimize_var_param.sh
# Instead of separate output and error log files, send the error output to the regular output log file
#$ -j y
# Where to put output/error files
#$ -o job_output_files
#$ -e job_output_files

echo "########################################################"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started date: "`date`
echo "##########################################################"

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript scripts/03_optimize_var_param.R --task_id $SGE_TASK_ID --n_tasks $SGE_TASK_LAST 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"