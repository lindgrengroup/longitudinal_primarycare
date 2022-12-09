# Author: Samvida S. Venkatesh
# Date: 10/11/2022

# Compare soft-clustering results from the two initialisations (true medoid, Euclidean distance-based medoid)

library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

# Read data ----

custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
names(custom_four_diverge) <- c("k1", "k2", "k3", "k4")

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

main_filepath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/"

soft_results <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    og_df <- read.table(paste0(main_filepath, p, "_", sx, 
                               "/soft_clustering_probs_", p, "_", sx, ".txt"),
                        sep = "\t", header = T, stringsAsFactors = F)
    colnames(og_df) <- c("eid", "k1_old", "k2_old", "k3_old", "k4_old")
    
    new_df <- read.table(paste0(main_filepath, p, "_", sx, 
                                "/medoid_init_soft_clustering_probs_", p, "_", sx, ".txt"),
                         sep = "\t", header = T, stringsAsFactors = F)
    colnames(new_df) <- c("eid", "k1_new", "k2_new", "k3_new", "k4_new")
    
    df <- full_join(og_df, new_df, by = "eid")
    return (df)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(soft_results) <- PHENOTYPES

# Plot old vs new soft-clust probabilities ----

oldVsNewPlot <- function (df, clustk) {
  sub_df <- df[, c("eid", paste0(clustk, "_old"), paste0(clustk, "_new"))]
  colnames(sub_df) <- c("eid", "old", "new")
  
  # Subset to individuals with > 0.75 or < 0.25 probabilities to distinguish clusters better
  sub_df <- sub_df %>%
    filter(old > 0.75 | old < 0.25 | new > 0.75 | new < 0.25)
  
  resplot <- ggplot(sub_df, aes(x = old, y = new)) +
    geom_point(colour = custom_four_diverge[clustk],
               fill = custom_four_diverge[clustk]) +
    labs(x = "old_dist_baselining", y = "new_dist_baselining",
         title = clustk)
  return (resplot)
}

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    per_clust <- lapply(c("k1", "k2", "k3", "k4"), function (clustk) {
      oldVsNewPlot(soft_results[[p]][[sx]], clustk)
    })
    arranged_plots <- ggarrange(plotlist = per_clust, nrow = 2, ncol = 2)
    png(paste0(main_filepath, p, "_", sx, "/plots/",
               p, "_", sx, "_compare_old_new_init.png"))
    print(arranged_plots)
    dev.off()
  })
})

