#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 31/10/22

# LiftOver FinnGen sumstats from hg38 to hg19

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N liftover_finngen_sumstats
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load liftOver/20210519

cd /well/lindgren/samvida/hormones_infertility/data/infertility_finngen
touch liftover_hg38Tohg19.log

cat pheno_list.txt | while read PHENO;
do
	printf "** Phenotype: ${PHENO}\n" >> liftover_hg38Tohg19.log 

	gunzip finngen_R7_${PHENO}.gz

	# Create bed file with chrX, pos0, pos1, rsid
	awk -v FS="\t" -v OFS="\t" 'NR > 1 {print "chr"$1, $2-1, $2, $5}' finngen_R7_${PHENO} \
	> tmp_finngen_R7_${PHENO}_hg38.bed
	awk -F "\t" 'BEGIN { OFS = FS } { if ($4 == "") $4 = $1":"$3; else $4 = $4; print }' \
	tmp_finngen_R7_${PHENO}_hg38.bed > finngen_R7_${PHENO}_hg38.bed

	liftOver finngen_R7_${PHENO}_hg38.bed \
	/well/lindgren/samvida/Resources/hg38ToHg19.over.chain.gz \
	finngen_R7_${PHENO}_hg19.bed \
	finngen_R7_${PHENO}_unlifted.bed

	# Log the % of SNPs succesfully lifted over
	hg38=$(cat finngen_R7_${PHENO}_hg38.bed | wc -l)
	printf "\t # SNPs in hg38: $hg38 \n" >> liftover_hg38Tohg19.log
	hg19=$(cat finngen_R7_${PHENO}_hg19.bed | wc -l)
	printf "\t # SNPs lifted over to hg19: $hg19 \n" >> liftover_hg38Tohg19.log
	perc_lifted=$(echo "scale=3;${hg19}/${hg38}" | bc)
	printf "\t percent successfully lifted: $(echo "scale=3;${hg19}/${hg38}*100" | bc) \n" >> liftover_hg38Tohg19.log

	# Create a new sumstats file by joining via rsid
	# rsid: 4th column in the lifted over bed file and 5th column in the original sumstats
	# Add header
	head -1 finngen_R7_${PHENO} > finngen_R7_${PHENO}_hg19_sumstats.txt
	join -1 4 -2 5 -t $'\t' -o 1.1,1.3,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13 \
	<(sort -k4 finngen_R7_${PHENO}_hg19.bed) <(sort -k5 finngen_R7_${PHENO}) >> finngen_R7_${PHENO}_hg19_sumstats.txt
	
	gzip finngen_R7_${PHENO}_hg19_sumstats.txt
	gzip finngen_R7_${PHENO}

	rm tmp_finngen_R7_${PHENO}_hg38.bed
	mv *.bed tmp_files	
done 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

