# Author: Samvida S. Venkatesh
# Date: 24/03/23

library(tidyverse)
theme_set(theme_bw())

# Read data ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/post_GWAS/lead_snps"

STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")

CLUSTERS <- c("k1", "k1_k2", "k1_k2_k3",
              "k1_no_FU_adjustment", "k1_k2_no_FU_adjustment", "k1_k2_k3_no_FU_adjustment")

dat <- lapply(STRATA, function (st) {
  res_list <- lapply(CLUSTERS, function (k) {
    res <- read.table(paste0(mainpath, "/all_strata_lead_snp_assocns/", st, "_", k,
                             "_lead_snp_assocns.txt"),
                      sep = "\t", header = T, stringsAsFactors = F)
    res$cluster <- k
    return (res)
  })
  res_list <- bind_rows(res_list)
  res_list$strata <- st
  return (res_list)
})
full_dat <- bind_rows(dat)

write.table(full_dat,
            paste0(mainpath, "/all_strata_lead_snp_assocns/all_strata_lead_snp_assocns.txt"),
            sep = "\t", quote = F, row.names = F)

# Wrangle for table ----

to_write <- full_dat %>%
  mutate(fuadj = ifelse(grepl("_no_FU_adjustment", cluster), "no", "adj"),
         cluster = gsub("_no_FU_adjustment", "", cluster),
         OR = exp(BETA),
         LCI = exp(BETA - 1.96*SE),
         UCI = exp(BETA + 1.96*SE),
         OR_write = paste0(signif(OR, 3), " (", signif(LCI, 3), " - ", signif(UCI, 3), ")"),
         CHRPOS_write = paste0(CHR, ":", POS)) %>%
  pivot_wider(names_from = fuadj,
              values_from = c(BETA, SE, PVALUE, OR, LCI, UCI, OR_write))

to_write <- to_write %>%
  filter(PVALUE_adj <= 5E-08 | PVALUE_no <= 5E-08) 

to_write <- to_write[, c("SNP", "CHRPOS_write", "AF_Tested",
                         "cluster", "strata", 
                         "OR_write_adj", "PVALUE_adj",
                         "OR_write_no", "PVALUE_no")]

write.table(to_write, "formatted_for_table.txt",
            sep = "\t", row.names = F, quote = F)

