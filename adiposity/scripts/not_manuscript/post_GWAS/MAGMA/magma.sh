#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 05/12/22
# ADAPTED FROM: Frederik Heymann Lassen

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N submit_magma
#$ -j y
#$ -V

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load gcccuda/2020b

STRATA="Weight_sex_comb"
PARAMETER="k1"
SAMPLE_SIZE=176408

readonly in_dir="/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/BOLT_results"
readonly out_dir="/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/post_GWAS/${STRATA}/MAGMA"

readonly sumstat="${in_dir}/${STRATA}_${PARAMETER}_final.txt"
readonly out_prefix="${out_dir}/${STRATA}_${PARAMETER}"
readonly snp_loc="${out_prefix}.snp_loc"

readonly magma_dir="/well/lindgren/samvida/Resources/MAGMA"
#readonly prefix_ref="${magma_dir}/auxiliary_files/reference/g1000_eur"
readonly ref_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/prs/hapmap/ld/unrel_kin_eur_10k"
readonly prefix_ref="${ref_dir}/short_merged_ukb_hapmap_rand_10k_eur"
readonly genes="${magma_dir}/auxiliary_files/genes/GRCh37/NCBI37.gene.loc"
readonly dbsnp="${magma_dir}/auxiliary_files/dbsnp/dbsnp151.synonyms"
readonly magma="${magma_dir}/./magma"

readonly rscript="magma.R"

# * Note: Only SNPs in reference panel and target are used for analysis and
# we should therefore generate a new panel using all imputed SNPs with info > 0.8
# * Note: MAGMA can not read in gzipped data

mkdir -p ${out_dir}

# Format gene gene loc file
cat ${sumstat} | awk '{print $1"\t"$2"\t"$3}' > ${snp_loc}

# Annotation Step (SNP to gene mapping)
if [ ! -f "${out_prefix}.genes.annot" ]; then
  set -x
  ${magma} \
    --annotate \
    --snp-loc "${snp_loc}" \
    --gene-loc "${genes}" \
    --out "${out_prefix}" 
  set +x
else
  >&2 echo "${out_prefix}.genes.annot already exist. Skipping.."
fi

# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
set -x
${magma} \
  --bfile "${prefix_ref}" \
  --gene-annot "${out_prefix}.genes.annot" \
  --pval "${sumstat}" "pval=9" "snp-id=1" "N=${SAMPLE_SIZE}" \
  --out "${out_prefix}"
set +x

module load R/4.0.5-foss-2020b

set -x
Rscript ${rscript} \
  --in_file "${out_prefix}.genes.out" \
  --out_file "${out_prefix}.txt" \
  --out_sep "\t"
set +x


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
