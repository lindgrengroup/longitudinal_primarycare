# Author: Samvida S. Venkatesh
# Date: 07/10/21

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(071021)

# Read data ----

covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/data/covariates.rds")
PHENO <- names(covars)
SEX_STRATA <- c("F", "M", "sex_comb")

rand_effs <- lapply(PHENO, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/random_effect_terms_",
                 p, ".rds"))
})
names(rand_effs) <- PHENO

# Plot distributions of each random effect term ----

coeff_distrib <- lapply(PHENO, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    res <- rand_effs[[p]][[sx]]
    res$sex <- sx
    return (res)
  })
  all_terms <- bind_rows(per_sex) %>% 
    pivot_longer(cols = -any_of(c("eid", "sex")),
                 names_to = "term", 
                 values_to = "coef_value")
  res_plot <- ggplot(all_terms, aes(x = sex, y = coef_value)) +
    facet_wrap(~term, scales = "free") +
    geom_violin(aes(fill = sex), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_y_continuous() +
    labs(x = "Sex strata", y = "Term coefficient", 
         title = paste0("Random effect distributions in: ", p)) +
    theme(legend.position = "none")
  return (res_plot)
})

pdf("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/random_effect_plots.pdf")
print(coeff_distrib)
dev.off()

# Test association of random effect coefficients with various covariates ----

covars_to_test <- c("baseline_age", "baseline_trait", "FU_n", "FUyrs")
assoc_tests <- lapply(PHENO, function (p) { 
  log_file <- paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/log_files/random_effect_associations_",
                     p, ".txt")
  per_sex <- lapply(SEX_STRATA, function (sx) {
    sink(log_file, append = T)
    cat(paste0("Running association tests in sex strata: ", sx, "\n"))
    sink()
    
    reff <- rand_effs[[p]][[sx]]
    colnames(reff)[grep("(Intercept)", colnames(reff))] <- "intrcpt"
    terms_to_test <- colnames(reff)[-which(colnames(reff) == "eid")]
    
    covar_dat <- covars[[p]]
    all_terms <- merge(reff, covar_dat, by = "eid")
    
    mod_results <- lapply(terms_to_test, function (tm) {
      all_covs <- lapply(covars_to_test, function (cov) {
          test_form <- formula(paste0(tm, " ~ ", cov))
          mod_res <- lm(test_form, data = all_terms)
          beta <- mod_res$coefficients[2]
          se <- sqrt(diag(vcov(mod_res)))[2]
          res <- data.frame(beta = beta, se = se,
                            term = tm, covariate = cov)
          return (res)
      })
      all_covs <- bind_rows(all_covs)
      return (all_covs)
    }) 
    mod_results <- bind_rows(mod_results)
    mod_results$sex_strata <- sx
    mod_results$phenotype <- p
    return (mod_results)
  })
  per_sex <- bind_rows(per_sex)
  write.table(per_sex, log_file, quote = F, sep = "\t", 
              row.names = F)
  return ()
})

# Plot association of random effect coefficients with covariates ----

assoc_plots <- lapply(PHENO, function (p) {
  
  per_sex <- lapply(SEX_STRATA, function (sx) {
    reff <- rand_effs[[p]][[sx]]
    reff_long <- reff %>% pivot_longer(cols = -any_of(c("eid")),
                                       names_to = "term", 
                                       values_to = "coef_value")
    covar_dat <- covars[[p]]
    all_terms <- merge(reff_long, covar_dat, by = "eid")
    
    bl_age_plot <- ggplot(all_terms, 
                          aes(x = baseline_age, y = coef_value)) +
      facet_wrap(~term, scales = "free") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Baseline age (years)", y = "Term coefficient", 
           title = paste0("Sex strata: ", sx, ", phenotype: ", p))
    
    bl_trait_plot <- ggplot(all_terms, 
                            aes(x = baseline_trait, y = coef_value)) +
      facet_wrap(~term, scales = "free") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Baseline trait", y = "Term coefficient", 
           title = paste0("Sex strata: ", sx, ", phenotype: ", p))
    
    fu_n_plot <- ggplot(all_terms, 
                        aes(x = FU_n, y = coef_value)) +
      facet_wrap(~term, scales = "free") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Number of follow-up measures", y = "Term coefficient", 
           title = paste0("Sex strata: ", sx, ", phenotype: ", p))
    
    fu_yrs_plot <- ggplot(all_terms, 
                          aes(x = FUyrs, y = coef_value)) +
      facet_wrap(~term, scales = "free") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Length of follow-up (years)", y = "Term coefficient", 
           title = paste0("Sex strata: ", sx, ", phenotype: ", p))
    
    return (list(bl_age_plot, bl_trait_plot, fu_n_plot, fu_yrs_plot))
  })
  
  pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/random_effect_associations_", 
             p, ".pdf"))
  print(per_sex)
  dev.off()
  
  return ()
})


