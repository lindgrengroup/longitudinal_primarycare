#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 31/10/22

# LiftOver ELGH sumstats from hg38 to hg19

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N liftover_elgh_sumstats
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load liftOver/20210519

cd /well/lindgren/samvida/hormones_infertility/data/infertility_elgh
touch liftover_hg38Tohg19.log

cat pheno_list.txt | while read PHENO;
do
	printf "** Phenotype: ${PHENO}\n" >> liftover_hg38Tohg19.log 

	gunzip 2022_06_30_GNH_binarytrait_GWAS_firth_ICD10__${PHENO}_ICD10__${PHENO}.regenie.gz

	# Create bed file with chrX, pos0, pos1, rsid
	awk -v FS=" " -v OFS="\t" 'NR > 1 {print "chr"$1, $2-1, $2, $3}' 2022_06_30_GNH_binarytrait_GWAS_firth_ICD10__${PHENO}_ICD10__${PHENO}.regenie \
	> ${PHENO}_hg38.bed
	
	liftOver ${PHENO}_hg38.bed \
	/well/lindgren/samvida/Resources/hg38ToHg19.over.chain.gz \
	${PHENO}_hg19.bed \
	${PHENO}_unlifted.bed

	# Log the % of SNPs succesfully lifted over
	hg38=$(cat ${PHENO}_hg38.bed | wc -l)
	printf "\t # SNPs in hg38: $hg38 \n" >> liftover_hg38Tohg19.log
	hg19=$(cat ${PHENO}_hg19.bed | wc -l)
	printf "\t # SNPs lifted over to hg19: $hg19 \n" >> liftover_hg38Tohg19.log
	perc_lifted=$(echo "scale=3;${hg19}/${hg38}" | bc)
	printf "\t percent successfully lifted: $(echo "scale=3;${hg19}/${hg38}*100" | bc) \n" >> liftover_hg38Tohg19.log

	# Create a new sumstats file by joining via rsid
	# rsid: 4th column in the lifted over bed file and 3rd column in the original sumstats
	# Add header
	head -1 2022_06_30_GNH_binarytrait_GWAS_firth_ICD10__${PHENO}_ICD10__${PHENO}.regenie \
	> 2022_06_30_GNH_binarytrait_GWAS_firth_ICD10__${PHENO}_hg19_sumstats.txt

	join -1 4 -2 3 -o 1.1,1.3,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 \
	<(sort -k 4 ${PHENO}_hg19.bed) <(sort -k 3 2022_06_30_GNH_binarytrait_GWAS_firth_ICD10__${PHENO}_ICD10__${PHENO}.regenie) \
	>> 2022_06_30_GNH_binarytrait_GWAS_firth_ICD10__${PHENO}_hg19_sumstats.txt
	
	gzip 2022_06_30_GNH_binarytrait_GWAS_firth_ICD10__${PHENO}_hg19_sumstats.txt
	gzip 2022_06_30_GNH_binarytrait_GWAS_firth_ICD10__${PHENO}_ICD10__${PHENO}.regenie

	mv *.bed tmp_files	
done 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
