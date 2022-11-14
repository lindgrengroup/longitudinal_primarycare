#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J cluster_sampling_scheme

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
