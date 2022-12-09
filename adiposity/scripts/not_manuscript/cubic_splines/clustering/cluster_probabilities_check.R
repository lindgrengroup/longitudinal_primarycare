# Author: Samvida S. Venkatesh
# Date: 03/05/22

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(RColorBrewer)

# Read files ----

PHENO <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

clust_results <- lapply(PHENO, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/",
                            p, "_", sx, "_remapped_clusters.txt"),
                     sep = "\t", header = T, stringsAsFactors = F)
    df$eid <- as.character(df$eid)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(clust_results) <- PHENO

# Get covariates to see whether clustering probabilities
# are associated with something other than just the spline models
covars <- lapply(PHENO, function (p) {
  readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[p]]
})
names(covars) <- PHENO

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220131_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

age_at_diag_matrix <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt",
                                 sep = "\t", header = T, stringsAsFactors = F)
colnames(age_at_diag_matrix) <- gsub("^X", "", colnames(age_at_diag_matrix))
age_at_diag_matrix <- age_at_diag_matrix[, 1:3]
age_at_diag_matrix$eid <- as.character(age_at_diag_matrix$eid)

CAT_COVARS <- "sex"
QUANT_COVARS <- c("baseline_age", "year_of_birth", "age_at_death",
                  "age_at_first_record", "age_at_last_record")
COVARS_PLOT <- c(CAT_COVARS, QUANT_COVARS)

# Wrangle data ----

# Add in age at first and last GP record to covariates
all_covars <- lapply(PHENO, function (p) {
  df <- covars[[p]]
  df$eid <- as.character(df$eid)
  
  df <- merge(df, general_covars, by = "eid")
  df <- merge(df, age_at_diag_matrix, by = "eid")

  return (df)
})
names(all_covars) <- PHENO

# Plots for max probability per individual -----

# Histogram
maxProbHist <- function (p, sx) {
  max_plot <- clust_results[[p]][[sx]] %>% 
    mutate(max_prob = select(., starts_with("k")) %>% do.call(pmax, .))
  
  res_plot <- ggplot(max_plot, aes(x = max_prob)) +
    geom_histogram() +
    expand_limits(y = 0)
  
  return (res_plot)
}

# vs covariates
maxProbVsCovars <- function (p, sx) {
  max_plot <- clust_results[[p]][[sx]] %>% 
    mutate(max_prob = select(., starts_with("k")) %>% do.call(pmax, .)) %>%
    select(all_of(c("eid", "max_prob"))) %>%
    # Cut max prob into 0.1 unit bins (<0.5, <0.6, etc.)
    mutate(max_prob = as.factor(as.character(paste0("<", 
                                                    plyr::round_any(max_prob, 0.1, f = floor)))))
  
  for_plot <- left_join(max_plot, all_covars[[p]],
                        by = "eid")
  
  plot_list <- lapply(QUANT_COVARS, function (cov) {
    res <- ggplot(for_plot %>% filter(!is.na(!!as.symbol(cov))), 
                  aes(x = max_prob, y = !!as.symbol(cov))) +
      geom_violin(aes(fill = max_prob), 
                  position = position_dodge(1)) +
      scale_fill_brewer(palette = "Dark2") +
      labs(x = "max prob of k", y = cov) +
      theme(legend.position = "none")
    return (res)
  })
  
  arranged_plots <- ggarrange(plotlist = plot_list, nrow = 3, ncol = 2)
  
  return (arranged_plots)
}

# # Plots for probability within clusters (facet wrap cluster) -----
# 
# clustProbVsCovars <- function (p, sx) {
#   clust_probs <- clust_results[[p]][[sx]] %>% 
#     select(starts_with("k") | "eid") %>%
#     pivot_longer(-eid, names_to = "cluster", values_to = "prob")
#   
#   for_plot <- left_join(clust_probs, all_covars[[p]],
#                         by = "eid")
#   
#   lapply(COVARS_PLOT, function (cov) {
#     res <- ggplot(for_plot %>% filter(!is.na(!!as.symbol(cov))), 
#                   aes(x = !!as.symbol(cov), y = prob)) +
#       facet_wrap(~cluster, ncol = 2) +
#       geom_violin(aes(fill = !!as.symbol(cov)), 
#                   position = position_dodge(1)) +
#       scale_fill_brewer(palette = "Dark2") +
#       labs(x = cov, y = "probability") +
#       theme(legend.position = "none")
#     
#     png(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/probabilities/",
#                p, "_", sx, "_cluster_prob_vs_", cov, ".png"))
#     print(res)
#     dev.off()
#     
#     return ()
#   })
#   return ()
# }

# Apply plotting functions ----

lapply(PHENO, function (p) {
  lapply(SEX_STRATA, function (sx) {
    
    png(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/probabilities/",
               p, "_", sx, "_max_prob.png"))
    print(maxProbHist(p, sx))
    dev.off()
    
    png(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/probabilities/",
               p, "_", sx, "_max_prob_vs_covars.png"))
    print(maxProbVsCovars(p, sx))
    dev.off()
    
    # The following function prints to png directly
    # clustProbVsCovars(p, sx)

  })
})

