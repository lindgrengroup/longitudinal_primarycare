# Author: Samvida S. Venkatesh
# Date: 19/07/22

library(tidyverse)

filepath_main <- "C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/2204_models/GWAS/"

# Read data ----

STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")

all_lead_snps <- lapply(STRATA, function (s) {
  res <- read.table(paste0(filepath_main,
                           "lead_snps/", s, "_lmm_intercepts_final.lead_snps.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  return (res)
})
names(all_lead_snps) <- STRATA

novel_refined_snps <- lapply(STRATA, function (s) {
  res <- read.table(paste0(filepath_main, s, "/lmm_intercepts/annotated_results_refined_novel_snps.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)
  return (res)
})
names(novel_refined_snps) <- STRATA

# Get percent of lead SNPs that are novel or refined ----

unique_lead_snps <- lapply(STRATA, function (s) {
  res <- all_lead_snps[[s]] %>% 
    mutate(chrpos = paste0(chromosome, ":", position)) %>%
    select(all_of(c("rsid", "chrpos"))) 
  return (res)
})
unique_lead_snps <- bind_rows(unique_lead_snps) %>%
  distinct()

unique_novel_refined_snps <- lapply(STRATA, function (s) {
  res <- novel_refined_snps[[s]] %>% 
    rename(rsid = SNP_og, chrpos = CHRPOS_og) %>%
    select(all_of(c("rsid", "chrpos", "status")))
  return (res)
})
unique_novel_refined_snps <- bind_rows(novel_refined_snps) %>%
  distinct()

# Calculate proportion of variance explained by novel or refined SNPs ----

prop_var_novel_refined <- lapply(STRATA, function (s) {
  df <- all_lead_snps[[s]]
  df <- df %>%
    mutate(varexp = (beta^2)*2*maf*(1 - maf),
           status = ifelse(rsid %in% novel_refined_snps[[s]]$SNP_og, 
                           "novel_refined", "reported"))
  vexp <- df %>% group_by(status) %>%
    summarise(varexp = paste0(round(sum(varexp)*100, 4), "%"), 
              n = n()) %>%
    mutate(strata = s)
  return (vexp)
})
prop_var_novel_refined <- bind_rows(prop_var_novel_refined)

# Grab sumstats for novel SNPs in each strata ----

novel_snps <- lapply(STRATA, function (s) {
  df <- all_lead_snps[[s]]
  novel_snps <- novel_refined_snps[[s]] %>% filter(status == "novel")
  novel_snps <- novel_snps$SNP_og
  res <- df %>% filter(rsid %in% novel_snps)
  res <- res %>% 
    mutate(strata = s,
           chrpos = paste0(chromosome, ":", position, ":", allele1, ":", allele2),
           beta_se = paste0(signif(beta, 3), " (", signif(se, 3), ")"))
  return (res[, c("strata", "rsid", "chrpos", "maf", "beta_se", "p")])
})
novel_snps <- bind_rows(novel_snps)

