# Author: Samvida S. Venkatesh
# Date: 27/07/21

library(mclust)
library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

set.seed(270721)

# Read files ----

# Parse in phenotype argument
args <- commandArgs(trailingOnly = T)
PHENO <- args[1]
SEX_STRATA <- c("F", "M", "sex_comb")

# Raw data
dat <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/adj_traits/adj_",
                      PHENO, ".rds"))

# Create clustering log file
log_file <- paste0(PHENO, "/baseline_clustering_log_", PHENO, ".txt")

# Get baseline adjusted trait value for clustering ----

cdat <- lapply(dat, function (df) {
  res <- df %>% group_by(eid) %>%
    arrange(eid, age_event) %>%
    summarise(baseline_adj_trait = first(adj_value))
  res <- res[, c("eid", "baseline_adj_trait")]
  return (res)
})

# Step through different values of clusters ----

MAXCLUST <- 6

mBIC <- lapply(SEX_STRATA, function (sx) {
  # Subset data for clustering
  for_clust <- cdat[[sx]]
  for_clust <- for_clust[, -which(colnames(for_clust) == "eid")]
  
  # Specify number of clusters
  all_BIC <- mclustBIC(for_clust, G = 1:MAXCLUST)
  # Pick best number of clusters
  mod_BIC <- Mclust(for_clust, x = all_BIC)
  # Save summary of chosen clusters
  sink(log_file, append = T)
  cat("Clustering in sex strata: ", sx, "\n")
  print(summary(mod_BIC, parameters = T))
  sink()
  return (list(model = mod_BIC, all_BIC = all_BIC))
})
names(mBIC) <- SEX_STRATA

# Visualise BIC
pdf(paste0(PHENO, "/baseline_mclust_plots_", PHENO, ".pdf"), onefile = T)
lapply(mBIC, function (x) plot(x$all_BIC))
dev.off()

# Assign individuals to clusters ---- 

clustered_dat <- lapply(SEX_STRATA, function (sx) {
  cluster_ids <- mBIC[[sx]]$model$classification
  res <- bind_cols(cdat[[sx]],
                   cluster = cluster_ids)
  write.table(res, 
              paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/splines/new_models/",
                     PHENO, "/baseline_clusters_", PHENO, "_", sx, ".txt"),
              sep = "\t", row.names = F, quote = F)
  return (res)
})
names(clustered_dat) <- SEX_STRATA

# Plot mean raw data in each group ---- 

clust_plots <- lapply(SEX_STRATA, function (sx) {
  # Add cluster number to raw data
  for_plot <- dat[[sx]][, c("eid", "age_event", "value")]
  for_plot <- subset(for_plot, for_plot$eid %in% clustered_dat[[sx]]$eid)
  for_plot$cluster <- clustered_dat[[sx]]$cluster[match(for_plot$eid,
                                                        clustered_dat[[sx]]$eid)]
  # Bin into age groups and calculate mean, S.E. of mean
  for_plot$age_bin <- cut(for_plot$age_event, seq(20, 80, by = 5))
  for_plot <- for_plot %>% 
    group_by(cluster, age_bin) %>% 
    summarise(mean = mean(value),
              count = n(),
              lwr = mean(value) - (1.96*sd(value)/sqrt(count)),
              upr = mean(value) + (1.96*sd(value)/sqrt(count)))
  for_plot$cluster <- as.factor(for_plot$cluster)
  
  kplot <- ggplot(for_plot, aes(x = age_bin, y = mean,
                                colour = cluster, fill = cluster,
                                group = cluster)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Age bin (years)", y = "Mean (95% CI of mean)",
         title = paste0("Sex strata: ", sx))
  
  return (kplot)
})
names(clust_plots) <- SEX_STRATA

pdf(paste0(PHENO, "/baseline_clustplot_", PHENO, ".pdf"), onefile = T)
print(clust_plots)
dev.off()

