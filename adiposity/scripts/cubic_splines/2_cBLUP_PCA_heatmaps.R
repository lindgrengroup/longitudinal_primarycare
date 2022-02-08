# Author: Samvida S. Venkatesh
# Date: 15/11/2021

# Submit with 5 cores

library(tidyverse)
library(pheatmap)
library(fastcluster)
theme_set(theme_bw())

set.seed(151121)

# Read files ----

args <- commandArgs(trailingOnly = T)
STRATA <- args[1]

p <- gsub("_.*", "", STRATA)
sx <- gsub(paste0(p, "_"), "", STRATA)

blups <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_blups_",
                        p, ".rds"))[[sx]]

# Calculate PCs and print scree plot ----

for_pca <- as.matrix(blups)
# Calculate PCs
pca_res <- prcomp(for_pca, scale = F)

# Plot scree plot
var_expl <- pca_res$sdev^2 / sum(pca_res$sdev^2)
var_expl <- data.frame(PC = paste0("PC", 1:length(var_expl)),
                       var = var_expl)
pca_plot <- ggplot(var_expl, aes(x = PC, y = var)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  labs(x = "PCs", y = "Proportion of variance explained",
       title = paste0("PCA of cubic spline BLUPs, phenotype: ", p,
                      " strata ", sx))

pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/PCA/blup_pca_screeplot_",
           p, "_", sx, ".pdf"))
print(pca_plot)
dev.off()

# Heatmaps of BLUPs ----

# RINT each coefficient for plotting
rint_x <- function (x) {
  return (qnorm((rank(x) - 0.5) / sum(!is.na(x))))
}

for_plot <- as.matrix(blups)
# Replace parentheses text with blanks
colnames(for_plot) <- gsub("\\s*\\([^\\)]+\\)", "", 
                           colnames(for_plot))

# Cluster with hclust from fastcluster 
clust_results <- hclust(dist(for_plot))

# For plotting, scale each column (RINT) to see differences better
for_plot <- apply(for_plot, 2, FUN = rint_x)

pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/PCA/blup_heatmaps_",
           p, "_", sx, ".pdf"))
pheatmap(for_plot, 
         cluster_rows = clust_results, 
         cluster_cols = F,
         scale = "none", 
         treeheight_col = 0, treeheight_row = 50,
         show_rownames = F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         legend = T, fontsize = 10,
         main = paste0("Cubic spline BLUPs in ", p, "_", sx))
dev.off()

