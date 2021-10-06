#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 14/05/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -t 1-12 -tc 4
#$ -N BMI_splines
#$ -j y

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" strata_filenames.txt`

module load R-bundle-Bioconductor/3.9-foss-2019a-R-3.6.0

Rscript 5_raw_slopes.R BMI $STRATA_NAME
