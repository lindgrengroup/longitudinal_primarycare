# Author: Samvida S. Venkatesh
# Date: 01/11/21

library(lme4)
library(tidyverse)
theme_set(theme_bw())

set.seed(011121)

# Read files ----

lmm_mods_path <- "" # REDACTED

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

slope_models <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0(lmm_mods_path, "/", p, "_full_model.rds"))
})
names(slope_models) <- PHENOTYPES

blups <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0(lmm_mods_path, "/", p, "_", sx, "_blups_full_model.txt"),
                      sep = "\t", header = T, stringsAsFactors = F)
    colnames(res) <- c("eid", "b0", "b1")
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(blups) <- PHENOTYPES

# QQ-plots ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    
    # QQ-plot for residuals
    for_resid_plot <- data.frame(resid = resid(slope_models[[p]][[sx]]))
    qq_resid <- ggplot(for_resid_plot, aes(sample = resid)) +
      stat_qq(size = 0.3) +
      stat_qq_line()
    ggsave(paste0(lmm_mods_path, "/plots/qq_residuals_", p, "_", sx, "_full_model.png"),
           qq_resid, height = 7, width = 7, units = "in")
    
    # QQ-plot for BLUPs
    qq_blup_b0 <- ggplot(blups[[p]][[sx]], aes(sample = b0)) +
      stat_qq(size = 0.3) +
      stat_qq_line()
    ggsave(paste0(lmm_mods_path, "/plots/qq_b0_BLUP_", p, "_", sx, "_full_model.png"),
           qq_blup_b0, height = 7, width = 7, units = "in")
    
    qq_blup_b1 <- ggplot(blups[[p]][[sx]], aes(sample = b1)) +
      stat_qq(size = 0.3) +
      stat_qq_line()
    ggsave(paste0(lmm_mods_path, "/plots/qq_b1_BLUP_", p, "_", sx, "_full_model.png"),
           qq_blup_b1, height = 7, width = 7, units = "in")
    
  })
})