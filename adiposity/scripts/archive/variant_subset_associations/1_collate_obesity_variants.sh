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

cd /well/lindgren/samvida/Resources

# Subset GWAS catalog associations with obesity-related traits
# $35 = MAPPED_TRAIT, $12=CHR, $13=POS, $22=SNP
awk -v FS="\t" -v OFS="\t" 'NR==FNR{a[$1]; next} {for (i in a) if ($35~i) print $12,$13,$13,$22}' \
GWASCatalog/gwascat_obesity_traits.txt GWASCatalog/all_associations_211102.txt \
> GWASCatalog/tmp_gwascat_obesity_associations.bed

# Get chr and pos for SNPs missing that information
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
R --file=./clean_gwascat_bedfile.R
rm GWASCatalog/tmp_gwascat_obesity_associations.bed

# Liftover hg38 SNPs to hg19 
module load liftOver/20210519

liftOver GWASCatalog/gwascat_obesity_associations_hg38.bed \
hg38ToHg19.over.chain.gz \
GWASCatalog/gwascat_obesity_associations_hg19.bed \
GWASCatalog/gwascat_obesity_associations_unlifted.bed

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

