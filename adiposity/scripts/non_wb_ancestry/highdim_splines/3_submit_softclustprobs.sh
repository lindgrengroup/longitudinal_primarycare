#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 3
#SBATCH -J replication_softclust

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/scripts/3_soft_cluster_probabilities.R \
--phenotype=${phenotype} \
--sex_strata=${ss} \
--nboots=100

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
