# Author: Samvida S. Venkatesh
# Date: 28/06/21

library(fpc)
library(tidyverse)
theme_set(theme_bw())

RANDOM_SEED <- 280621
set.seed(RANDOM_SEED)

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--sex_strata", required = TRUE,
                    help = "Sex strata")

p <- args$phenotype
sx <- args$sex_strata

plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/plots/stability/"
resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/logs/"

KTEST <- 2:10 # number of clusters to evaluate
NBOOTS <- 50 # number of times to run bootstrap clusters to evaluate assignment

NSAMPLES <- 1000 # number of individuals to sample each iteration
NITERS <- 9 # number of times to draw samples from population

# Read files ----

blups <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/spline_models/",
                                    p, "_", sx, "_blups_full_model.txt"), 
                             sep = "\t", header = T, stringsAsFactors = F)

# Use bootstrap to select number of clusters based on stability metrics ----

dat <- as.matrix(blups[, -1])
rownames(dat) <- blups$eid

all_iters <- lapply(1:NITERS, function (ni) {
  sub_dat <- dat[sample(1:nrow(dat), NSAMPLES, replace = F), ]
  sub_nclust <- nselectboot(sub_dat, 
                            B = NBOOTS, 
                            clustermethod = claraCBI, 
                            k = KTEST)
  stab_mat <- 1 - sub_nclust$stab
  to_plot <- data.frame(iter_index = ni,
                        nk = 1:ncol(stab_mat),
                        mean_stab = colMeans(stab_mat), 
                        sd_stab = apply(stab_mat, 2, function (x) sd(x)))
  return (to_plot)
})
all_iters <- bind_rows(all_iters)
saveRDS(all_iters, 
        paste0(resdir, "clustering_stability_", p, "_", sx, ".rds"))

res_plot <- ggplot(all_iters %>% filter(nk >= 2), 
                   aes(x = nk, y = mean_stab)) +
  facet_wrap(~iter_index, nrow = 3, ncol = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_stab - 1.96*sd_stab,
                    ymax = mean_stab + 1.96*sd_stab), width = 0.7) +
  geom_hline(yintercept = 0.6, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.85, linetype = "dashed", colour = "blue") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = KTEST, labels = KTEST) +
  labs(x = "# clusters", y = "stability index")

# Visualise stability
pdf(paste0(plotdir, "clustering_stability_", p, "_", sx, ".pdf"))
plot(res_plot)
dev.off()
