#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 4
#$ -N cluster_sampling_scheme
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/scripts/2_sample_clustering_scheme.R \
--phenotype=${phenotype} \
--ss=${ss} \
--K=${nclust} \
--L=${lmin} \
--M=${myrs}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
