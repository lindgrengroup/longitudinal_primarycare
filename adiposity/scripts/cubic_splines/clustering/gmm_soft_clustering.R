# Author: Samvida S. Venkatesh
# Date: 28/06/21

library(mclust)
library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

set.seed(280621)

# Read files ----

# # Parse in phenotype argument
# args <- commandArgs(trailingOnly = T)
# STRATA <- args[1]
# 
# # Get coefficients for clustering (spline full model BLUPs)
# cdat <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/",
#                           STRATA, "_blups_full_model.txt"), 
#                    sep = "\t", header = T, stringsAsFactors = F)
# 
# # Step through different values of clusters ----
# 
# MAXCLUST <- 10
# 
# for_clust <- cdat[, -which(colnames(cdat) == "eid")]
# for_clust <- for_clust %>% mutate(across(everything(), as.numeric))
# # Specify number of clusters
# all_BIC <- mclustBIC(for_clust, G = 1:MAXCLUST)
# 
# # Visualise BIC
# png(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/",
#            STRATA, "_BIC.png"))
# plot(all_BIC)
# dev.off()

#### STOP HERE TO PICK BEST NUMBER OF CLUSTERS #####

## Continue

library(argparse)
library(mclust)
library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

set.seed(280621)

# Read in arguments ----

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Strata for trait and sex, ex. Weight_F")
parser$add_argument("--nClust", required=TRUE,
                    help = "# clusters")
args <- parser$parse_args()

# Read data ----

# Get coefficients for clustering (spline full model BLUPs)
cdat <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/",
                          args$strata, "_blups_full_model.txt"), 
                   sep = "\t", header = T, stringsAsFactors = F)

# Create clustering log file
log_file <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/logs/",
                   args$strata, "_spline_full_models.txt")

# Create clustering model ----

for_clust <- cdat[, -which(colnames(cdat) == "eid")]
for_clust <- for_clust %>% mutate(across(everything(), as.numeric))
best_mod <- Mclust(for_clust, G = args$nClust)

# Save summary of chosen model
sink(log_file, append = T)
print(summary(best_mod, parameters = T))
sink()

# Add cluster probabilities and assigned membership to dataframe of ids 
probs <- data.frame(best_mod$z)
colnames(probs) <- paste0("k", 1:ncol(probs))

id_clusts <- data.frame(eid = cdat$eid)
id_clusts <- bind_cols(id_clusts, probs)

assignments <- data.frame(assigned_k = paste0("k", best_mod$classification))
id_clusts <- bind_cols(id_clusts, assignments)

write.table(id_clusts, 
            paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/",
                   args$strata, "_spline_full_models.txt"),
            sep = "\t", row.names = F, quote = F)

# Plots ----

# Plot probability of cluster membership 
kprobs <- id_clusts %>% 
  pivot_longer(cols = -all_of(c("eid", "assigned_k")), 
               names_to = "cluster", values_to = "probability")

kplot <- ggplot(kprobs, aes(x = probability)) +
  facet_wrap(~cluster, nrow = 2) +
  geom_histogram() 

png(filename = paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/",
                      args$strata, "_cluster_probabilities.png"))
print(kplot)
dev.off()

# Plot coefficient correlations coloured by cluster 

termCorrPlot <- function (mat, t1, t2) {
  # Plot
  mat$assigned_k <- factor(as.character(mat$assigned_k))
  res_plot <- ggplot(mat,
                     aes(x = !!as.symbol(t1), y = !!as.symbol(t2),
                         colour = assigned_k)) +
    geom_point(size = 0.3) +
    scale_color_brewer(palette = "Set1")
  return (res_plot)
}

plot_coeffs <- cdat
plot_coeffs$assigned_k <- id_clusts$assigned_k[match(plot_coeffs$eid,
                                                     id_clusts$eid)]

TNAMES <- colnames(plot_coeffs)[!colnames(plot_coeffs) %in% c("eid", "assigned_k")]
plot_list <- lapply(1:(length(TNAMES)-1), function (i) {
  res <- lapply((i+1):length(TNAMES), function (j) {
    png(filename = paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/",
                          args$strata, "_", "coeff_corrs_",
                          i, "_", j, "_cluster_assignment.png"))
    print(termCorrPlot(plot_coeffs, TNAMES[i], TNAMES[j]))
    dev.off()
  })
  return ()
})

