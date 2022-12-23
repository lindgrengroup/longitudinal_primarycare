# Author: Samvida S. Venkatesh
# Date: 07/10/21

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(071021)

# Read data ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1

blups <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm_blups_",
                 p, ".rds"))
})
names(blups) <- PHENOTYPES

SEX_STRATA <- c("F", "M")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

PCs <- paste0("PC", 1:21)
COVARS <- c("baseline_age", "age_sq", "data_provider")

# Plot distributions of BLUPs ----

blup_distrib <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    res <- blups[[p]][[sx]]
    res$eid <- rownames(res)
    res <- res %>% pivot_longer(cols = -eid, 
                                names_to = "term",
                                values_to = "BLUP")
    res_plot <- ggplot(res, aes(x = BLUP)) +
      facet_wrap(~term, scales = "free") +
      geom_density() +
      labs(x = "BLUP", title = paste0("BLUP distributions in ",
                                      sx, " phenotype: ", p)) +
      theme(legend.position = "none")
    return (res_plot)
  })
  pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/lmm_blup_distrib_",
             p, ".pdf"))
  print(per_sex)
  dev.off()
  return ()
})

# Plot association of fixed+random effects with covariates ----

assoc_plots <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    bp <- blups[[p]][[sx]]
    blup_terms <- colnames(bp)
    bp$eid <- rownames(bp)
    covar_dat <- covars[[p]]
    all_terms <- merge(bp, covar_dat, by = "eid")
    # Randomly select subset of ids to plot,
    # as there is otherwise too much data to plot
    nsample <- min(nrow(all_terms), 500)
    all_terms <- all_terms[sample(nrow(all_terms), nsample,
                                  replace = F), ]
    all_terms <- pivot_longer(all_terms, 
                              cols = all_of(blup_terms),
                              names_to = "term",
                              values_to = "BLUP")
    
    bl_age_plot <- ggplot(all_terms, 
                          aes(x = baseline_age, y = BLUP)) +
      facet_wrap(~term, scales = "free") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Baseline age (years)", y = "BLUP", 
           title = paste0("BLUP vs baseline age in ",
                          sx, " phenotype: ", p))
    
    bl_trait_plot <- ggplot(all_terms, 
                            aes(x = baseline_trait, y = BLUP)) +
      facet_wrap(~term, scales = "free") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Baseline trait", y = "BLUP", 
           title = paste0("BLUP vs baseline trait in ",
                          sx, " phenotype: ", p))
    
    fu_n_plot <- ggplot(all_terms, 
                        aes(x = FU_n, y = BLUP)) +
      facet_wrap(~term, scales = "free") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Number of follow-up measures", y = "BLUP", 
           title = paste0("BLUP vs # follow-ups in ",
                          sx, " phenotype: ", p))
    
    fu_yrs_plot <- ggplot(all_terms, 
                          aes(x = FUyrs, y = BLUP)) +
      facet_wrap(~term, scales = "free") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Length of follow-up (years)", y = "BLUP", 
           title = paste0("BLUP vs # years of follow-up in ",
                          sx, " phenotype: ", p))
    return (list(bl_age_plot, bl_trait_plot, fu_n_plot, fu_yrs_plot))
  }) 
  pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/lmm_blup_assocns_",
             p, ".pdf"))
  print(per_sex)
  dev.off()
  return ()
})


