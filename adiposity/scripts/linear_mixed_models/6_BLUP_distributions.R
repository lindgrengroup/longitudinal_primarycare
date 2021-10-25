# Author: Samvida S. Venkatesh
# Date: 07/10/21

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(071021)

# Read data ----

covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/covariates.rds")
PHENOTYPES <- names(covars)
SEX_STRATA <- c("F", "M", "sex_comb")

blups <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm_blups_all_strata.rds")

# Plot distributions of BLUPs ----

blup_distrib <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    res <- blups[[p]][[sx]]
    res$sex <- sx
    return (res)
  })
  per_sex <- bind_rows(per_sex)
  per_sex$pheno <- p
  return (per_sex)
})
all_terms <- bind_rows(blup_distrib) 
blup_distrib_plot <- ggplot(all_terms, aes(x = sex, y = lmm_slope)) +
  facet_wrap(~pheno, scales = "free") +
  geom_violin(aes(fill = sex), position = position_dodge(1)) +
  geom_boxplot(width = 0.1) + 
  scale_y_continuous() +
  labs(x = "Sex strata", y = "BLUP", 
       title = "BLUP distributions") +
  theme(legend.position = "none")

# Test association of BLUPs with various covariates ----

covars_to_test <- c("baseline_age", "baseline_trait", "FU_n", "FUyrs")
log_file <- "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/log_files/BLUP_associations.txt"
assoc_tests <- lapply(PHENOTYPES, function (p) { 
  per_sex <- lapply(SEX_STRATA, function (sx) {
    sink(log_file, append = T)
    cat(paste0("** Phenotype and Strata: ", p, ", ", sx, "\n"))
    sink()
    
    bp <- blups[[p]][[sx]]
    covar_dat <- covars[[p]]
    full_dat <- merge(bp, covar_dat, by = "eid")
    
    mod_results <- lapply(covars_to_test, function (cov) {
      test_form <- formula(paste0("lmm_slope ~ ", cov))
      mod_res <- lm(test_form, data = full_dat)
      beta <- mod_res$coefficients[2]
      se <- sqrt(diag(vcov(mod_res)))[2]
      res <- data.frame(beta = beta, se = se,
                        covariate = cov)
      return (res)
    }) 
    mod_results <- bind_rows(mod_results)
    mod_results$sex_strata <- sx
    mod_results$phenotype <- p
    return (mod_results)
  })
  per_sex <- bind_rows(per_sex)
  return (per_sex)
})
assoc_tests <- bind_rows(assoc_tests)
write.table(assoc_tests, log_file, quote = F, sep = "\t", 
            row.names = F)

# Plot association of random effect coefficients with covariates ----

assoc_plots <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    bp <- blups[[p]][[sx]]
    covar_dat <- covars[[p]]
    all_terms <- merge(bp, covar_dat, by = "eid")
    all_terms$sex <- sx
    return (all_terms)
  }) 
  to_plot <- bind_rows(per_sex)
  bl_age_plot <- ggplot(to_plot, 
                        aes(x = baseline_age, y = lmm_slope)) +
    facet_wrap(~sex, scales = "free") +
    geom_point() + 
    geom_smooth(method = "lm", formula = y~x) +
    labs(x = "Baseline age (years)", y = "BLUP", 
         title = paste0("Phenotype: ", p))
  
  bl_trait_plot <- ggplot(to_plot, 
                          aes(x = baseline_trait, y = lmm_slope)) +
    facet_wrap(~sex, scales = "free") +
    geom_point() + 
    geom_smooth(method = "lm", formula = y~x) +
    labs(x = "Baseline trait", y = "BLUP", 
         title = paste0("Phenotype: ", p))
  
  fu_n_plot <- ggplot(to_plot, 
                      aes(x = FU_n, y = lmm_slope)) +
    facet_wrap(~sex, scales = "free") +
    geom_point() + 
    geom_smooth(method = "lm", formula = y~x) +
    labs(x = "Number of follow-up measures", y = "BLUP", 
         title = paste0("Phenotype: ", p))
  
  fu_yrs_plot <- ggplot(to_plot, 
                        aes(x = FUyrs, y = lmm_slope)) +
    facet_wrap(~sex, scales = "free") +
    geom_point() + 
    geom_smooth(method = "lm", formula = y~x) +
    labs(x = "Length of follow-up (years)", y = "BLUP", 
         title = paste0("Phenotype: ", p))
  
  return (list(bl_age_plot, bl_trait_plot, fu_n_plot, fu_yrs_plot))
})


pdf("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/lmm_blup_plots.pdf")
print(blup_distrib_plot)
print(assoc_plots)
dev.off()
