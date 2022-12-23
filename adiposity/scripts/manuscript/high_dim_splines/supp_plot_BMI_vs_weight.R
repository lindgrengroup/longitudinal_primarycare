# Author: Samvida S. Venkatesh
# Date: 07/12/2022

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(071222)

# Read data ----

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
hidim_mods_path <- "" # REDACTED

SEX_STRATA <- c("F", "M", "sex_comb")

df_list <- lapply(SEX_STRATA, function (sx) {
  bmi_probs <- read.table(paste0(hidim_mods_path, "/clustering/BMI_", 
                                 sx, "/soft_clustering_probs_BMI_", sx, ".txt"),
                          sep = "\t", header = T, stringsAsFactors = F)
  bmi_probs <- bmi_probs %>%
    mutate(eid = as.character(eid),
           k1_k2 = k1 + k2,
           k1_k2_k3 = k1 + k2 + k3) %>%
    select(all_of(c("eid", "k1", "k1_k2", "k1_k2_k3")))
  colnames(bmi_probs) <- c("eid", "BMI_k1", "BMI_k1_k2", "BMI_k1_k2_k3")
  
  weight_probs <- read.table(paste0(hidim_mods_path, "/clustering/Weight_", 
                                 sx, "/soft_clustering_probs_Weight_", sx, ".txt"),
                          sep = "\t", header = T, stringsAsFactors = F)
  weight_probs <- weight_probs %>%
    mutate(eid = as.character(eid),
           k1_k2 = k1 + k2,
           k1_k2_k3 = k1 + k2 + k3) %>%
    select(all_of(c("eid", "k1", "k1_k2", "k1_k2_k3")))
  colnames(weight_probs) <- c("eid", "Weight_k1", "Weight_k1_k2", "Weight_k1_k2_k3")
  
  df <- inner_join(bmi_probs, weight_probs, by = "eid")
  return (df)
})
names(df_list) <- SEX_STRATA

# Calculate correlation between bmi and weight probs in each strata ----

lapply(SEX_STRATA, function (sx) {
  res <- df_list[[sx]]
  print(paste0("Sex strata: ", sx, " N individuals: ", nrow(res)))
  print(paste0("k1 corr: ", signif(cor(res$BMI_k1, res$Weight_k1), 3)))
  print(paste0("k1_k2 corr: ", signif(cor(res$BMI_k1_k2, res$Weight_k1_k2), 3)))
  print(paste0("k1_k2_k3 corr: ", signif(cor(res$BMI_k1_k2_k3, res$Weight_k1_k2_k3), 3)))
})

