# Author: Samvida S. Venkatesh
# Date: 07/12/2022

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(071222)

# Read data ----

lmm_mods_path <- "" # REDACTED

SEX_STRATA <- c("F", "M", "sex_comb")

df_list <- lapply(SEX_STRATA, function (sx) {
  bmi_blups <- read.table(paste0(lmm_mods_path, "/BMI_", sx, "_blups_full_model.txt"),
                          sep = "\t", header = T, stringsAsFactors = F)
  bmi_blups$eid <- as.character(bmi_blups$eid)
  colnames(bmi_blups) <- c("eid", "BMI_b0", "BMI_b1")
  
  weight_blups <- read.table(paste0(lmm_mods_path, "/Weight_", sx, "_blups_full_model.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
  weight_blups$eid <- as.character(weight_blups$eid)
  colnames(weight_blups) <- c("eid", "Weight_b0", "Weight_b1")
  
  df <- inner_join(bmi_blups, weight_blups, by = "eid")
  return (df)
})
names(df_list) <- SEX_STRATA

# Calculate correlation between bmi_b1 and weight_b1 in each strata ----

lapply(SEX_STRATA, function (sx) {
  print(paste0("Sex strata: ", sx, " N individuals: ", nrow(df_list[[sx]]), 
               " - b1 corr: ", 
               signif(cor(df_list[[sx]]$BMI_b1, df_list[[sx]]$Weight_b1), 3)))
})

