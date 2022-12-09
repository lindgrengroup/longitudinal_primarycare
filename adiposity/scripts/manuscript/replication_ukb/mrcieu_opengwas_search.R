# Author: Samvida S. Venkatesh
# Date: 23/08/22

# Perform look-ups of SNPs in published GWAS from MRCIEU GWAS catalog
library(ieugwasr)
library(tidyverse)

# Read SNP list 
snps_check <- read.table("snps_to_replicate.txt", sep = "\t", header = F,
                         stringsAsFactors = F)$V1

# List the non-UKB datasets to examine for replication datasets 
blist <- batches()
# Can't use the following because they are UKB data:
# ukb-a, ukb-b, ukb-d, ukb-e
# Can't use the following because they don't contain anthropometric measures:
# eqtl-a
# finn-b
# met-a, met-b, met-c, met-d
# prot-a, prot-b, prot-c
# ubm-a
bkeep <- c("bbj-a", "ebi-a", "ieu-a", "ieu-b")

# Perform phewas for batches with anthropometric measures
phewas_res <- phewas(variants = snps_check, pval = 1E-3, batch = bkeep)
phewas_res <- split(phewas_res, f = phewas_res$rsid)

phewas_res <- lapply(phewas_res, function (df) {
  return (df %>% arrange(p))
})

View(phewas_res[["rs2861761"]])

