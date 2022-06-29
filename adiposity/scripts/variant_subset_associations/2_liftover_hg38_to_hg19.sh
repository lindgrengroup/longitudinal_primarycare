#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/04/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -N liftover_hg38_to_hg19
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load liftOver/20210519

cd /well/lindgren/samvida/Resources

liftOver GWASCatalog/gwascat_obesity_associations.bed \
hg38ToHg19.over.chain.gz \
gwascat_obesity_associations_hg19.bed \
gwascat_obesity_associations_unlifted.bed

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
