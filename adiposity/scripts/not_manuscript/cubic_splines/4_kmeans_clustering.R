# Author: Samvida S. Venkatesh
# Date: 01/02/2022

library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

set.seed(010222)

# Read files ----

PHENOTYPES <- c("BMI", "WC", "Weight", "WHR")
SEX_STRATA <- c("F", "M", "sex_comb")

blups <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_blups_",
                 p, ".rds"))
})
names(blups) <- PHENOTYPES

covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]
CLUST_COVARS <- c("baseline_age", "age_sq", "year_of_birth",
                  "baseline_trait", "FUyrs", "FU_n")

# Wrangle data to add covariates to clustering ----

to_cluster <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) { 
    coeffs <- blups[[p]][[sx]]
    coeffs$eid <- rownames(coeffs)
    relevant_covars <- covars[[p]]
    df <- merge(coeffs, relevant_covars[, c("eid", CLUST_COVARS)],
                 by = "eid")
    res <- as.matrix(df[, -1])
    rownames(res) <- df$eid
    return (res)
  }) 
  names(per_sex) <- SEX_STRATA
  return (per_sex)
})
names(to_cluster) <- PHENOTYPES

# Step through different values of clusters ----

MAXCLUST <- 20

all_wss <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    
    kmods <- lapply(1:MAXCLUST, function (ki) {
      res <- kmeans(to_cluster[[p]][[sx]], centers = ki, iter.max = 1000)
      wss <- res$tot.withinss
      return (wss)
    })
    
    dat_plot <- data.frame(k = 1:MAXCLUST,
                           wss = unlist(kmods))
    to_plot <- ggplot(dat_plot, aes(x = k, y = wss)) +
      geom_point() +
      geom_line() +
      scale_x_continuous(breaks = 1:MAXCLUST) +
      labs(x = "Number of clusters k", y = "Total within sum of squares",
           title = paste0("k-means clustering in pheno: ", p, " strata: ", sx))
    # Visualise BIC
    pdf(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/plots/cubic_splines/clustering/kmeans_wss_",
               p, "_", sx, ".pdf"))
    print(to_plot)
    dev.off()
    return ()
  })
  return ()
})

###### BREAK HERE TO DETERMINE THE NUMBER OF CLUSTERS FROM WSS PLOTS ######
## 5 seems reasonable

# Pick best number of clusters
mod_kmeans <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    strata_nm <- paste0(p, "_", sx)
    nclust <- 5
    res <- kmeans(to_cluster[[p]][[sx]], 
                  centers = nclust,
                  iter.max = 10000,
                  algorithm = "Lloyd")
    
    # Save summary of chosen clusters
    sink(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/log_files/kmeans_cubic_splines_",
                p, ".txt"), append = T)
    cat("Clustering in sex strata: ", sx, "\n")
    print(res$centers)
    sink()
    
    to_save <- data.frame(eid = rownames(blups[[p]][[sx]]),
                          k = res$cluster)
    write.table(to_save, 
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/results/cubic_splines/clustering_with_covars_",
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    return ()
  })
  return ()
})

