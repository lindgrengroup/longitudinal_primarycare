# Author: Samvida S. Venkatesh
# Date: 28/06/21

library(mclust)
library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

set.seed(280621)

# Read files ----

PHENOTYPES <- c("BMI", "WC", "Weight", "WHR")
SEX_STRATA <- c("F", "M", "sex_comb")

blups <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/adj_FUyrs_not_genetic_PCs/adj_FUyrs_cubic_spline_blups_",
                 p, ".rds"))
})
names(blups) <- PHENOTYPES

NCLUSTERS_TABLE <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/log_files/visually_determined_nclusters.txt",
                              sep = "\t", header = T, stringsAsFactors = F) 

# Step through different values of clusters ----

MAXCLUST <- 9

all_BIC <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    res <- mclustBIC(blups[[p]][[sx]], G = 1:MAXCLUST)
    # Visualise BIC
    pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/clustering/clustering_BIC_",
               p, "_", sx, ".pdf"))
    plot(res)
    dev.off()
    return ()
  })
  return ()
})

###### BREAK HERE TO DETERMINE THE NUMBER OF CLUSTERS FROM BIC PLOTS ######

# Pick best number of clusters
mod_BIC <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    strata_nm <- paste0(p, "_", sx)
    res <- Mclust(blups[[p]][[sx]], 
                  G = NCLUSTERS_TABLE$nclusters[NCLUSTERS_TABLE$strata == strata_nm])
    
    # Save summary of chosen clusters
    sink(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/log_files/cubic_spline_clustering_",
                p, ".txt"), append = T)
    cat("Clustering in sex strata: ", sx, "\n")
    print(summary(res, parameters = T))
    sink()
    
    to_save <- data.frame(eid = rownames(blups[[p]][[sx]]),
                          k = res$classification)
    write.table(to_save, 
                paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/adj_FUyrs_not_genetic_PCs/adj_FUyrs_clusters_",
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    return ()
  })
  return ()
})

