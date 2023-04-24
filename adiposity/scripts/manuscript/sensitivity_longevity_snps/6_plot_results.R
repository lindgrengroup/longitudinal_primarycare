# Author: Samvida S. Venkatesh
# Date: 11/04/2023

library(tidyverse)
theme_set(theme_bw())

# Read data ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/grepped_sumstats"

PHENO <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")
PARAMETERS <- c("b1", "k1", "k1_k2", "k1_k2_k3")

assoc_dat <- lapply(PHENO, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    per_par <- lapply(PARAMETERS, function (pr) {
      res <- read.table(paste0(mainpath, "/", p, "_", sx, "_", pr, "_longevity_rsids.txt"), sep = "\t", header = T, 
                        comment.char = "@", stringsAsFactors = F)
      res$pheno <- p
      res$sex_strata <- sx
      res$parameter <- pr
      return (res)
    })
    per_par <- bind_rows(per_par)
    return (per_par)
  })
  res_list <- bind_rows(res_list)
  return (res_list)
})
assoc_dat <- bind_rows(assoc_dat)

# Examine results ----

ntests <- length(unique(assoc_dat$SNP))
PTHRESH <- 0.05/ntests
GWS_THRESH <- 5E-08

write.table(assoc_dat, paste0(mainpath, "/all_results_grepped_sumstats.txt"),
            sep = "\t", row.names = F, quote = F)

# Pivot wider to write table ----

to_write <- assoc_dat %>% 
  mutate(MAF = ifelse(AF_Tested < 0.5, AF_Tested, 1 - AF_Tested),
         strata = paste0(pheno, "_", sex_strata)) %>%
  pivot_wider(names_from = parameter,
              values_from = c(BETA, SE, PVALUE)) %>%
  arrange(CHR, POS, strata)

to_write <- to_write[, c("SNP", "CHR", "POS", "MAF", "strata",
                         "BETA_b1", "SE_b1", "PVALUE_b1",
                         "BETA_k1", "SE_k1", "PVALUE_k1",
                         "BETA_k1_k2", "SE_k1_k2", "PVALUE_k1_k2",
                         "BETA_k1_k2_k3", "SE_k1_k2_k3", "PVALUE_k1_k2_k3")]
write.table(to_write, "wide_for_print.txt",
            sep = "\t", row.names = F, quote = F)
