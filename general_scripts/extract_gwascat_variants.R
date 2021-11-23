# Author: Samvida S. Venkatesh
# Date: 01/11/2021

library(tidyverse)

# Read list of Pubmed IDs ----

me_gwas <- read.table("largest_gwas_metabolic_endocrine.txt",
                            sep = "\t", header = T, quote = "")
pub_gwas <- me_gwas[!is.na(me_gwas$Pubmed.ID), ]

all_assocns <- read.table("all_associations_211102.txt",
                          sep = "\t", header = T, quote = "")

# Get index variants from pubmed id ----

tmp_assocns <- all_assocns %>% 
  mutate(match_col = paste0(PUBMEDID, DISEASE.TRAIT, sep = "_"))
tmp_gwas <- pub_gwas %>% 
  mutate(match_col = paste0(Pubmed.ID, Trait.or.disorder, sep = "_"))

me_assocns <- tmp_assocns %>%
  filter(match_col %in% tmp_gwas$match_col) %>%
  select(PUBMEDID, FIRST.AUTHOR, DISEASE.TRAIT, REGION,
         CHR_ID, CHR_POS, STRONGEST.SNP.RISK.ALLELE,
         SNPS, MERGED, SNP_ID_CURRENT, RISK.ALLELE.FREQUENCY,
         P.VALUE, PVALUE_MLOG, P.VALUE..TEXT., 
         OR.or.BETA, X95..CI..TEXT.)

write.table(me_assocns, "metabo_endo_SNPs_largest_gwas.txt",
            sep = "\t", row.names = F, quote = F)
