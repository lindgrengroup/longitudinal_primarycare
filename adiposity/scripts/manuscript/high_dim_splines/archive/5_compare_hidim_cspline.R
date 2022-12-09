# Author: Samvida S. Venkatesh
# Date: 26/05/2022

library(argparse)
library(tidyverse)
library(zoo)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--sex_strata", required = TRUE,
                    help = "Sex strata")

args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$sex_strata

plotdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/plots/final_clusters/"
resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/clustering/"

# Load data ----

## Clustering results

hidim_clusts <- read.table(paste0(resdir, "assigned_clusters_", PHENO, "_", SEX_STRATA, 
                                  ".txt"),
                           sep = "\t", header = T, stringsAsFactors = F)
hidim_clusts <- hidim_clusts %>% 
  rename(hidim_k = clust) %>%
  mutate(eid = as.character(eid),
         hidim_k = as.character(hidim_k))

cspline_clusts <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/",
                                    PHENO, "_", SEX_STRATA, "_spline_full_models.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
cspline_clusts <- cspline_clusts %>% 
  select(all_of(c("eid", "assigned_k"))) %>%
  mutate(eid = as.character(eid),
         assigned_k = gsub("k", "", assigned_k))

NCLUST <- length(unique(hidim_clusts$hidim_k)) 

## Raw data and metadata

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[[PHENO]]
# Covariates
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]
covars$eid <- as.character(covars$eid)

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Remap clusters for csplines to line up somewhat with hidim splines ----

remap_table <- data.frame(old_k = as.character(c(2, 1, 5, 3, 4)),
                          new_k = as.character(c(1:5)))
cspline_clusts$cspline_k <- 
  remap_table$new_k[match(cspline_clusts$assigned_k,
                          remap_table$old_k)]

# Create cross-sectional quintile assignments ----

quintile_dat <- dat %>% group_by(eid) %>%
  arrange(event_dt, .by_group = T) %>%
  summarise(baseline_value = first(value),
            mean_value = mean(value)) %>%
  mutate(baseline_k = as.character(ntile(baseline_value, NCLUST)),
         mean_k = as.character(ntile(mean_value, NCLUST))) %>%
  select(all_of(c("eid", "baseline_k", "mean_k")))

# Make 1 the highest quintile and 5 the lowest
remap_table <- data.frame(old_k = as.character(c(1:5)),
                          new_k = as.character(c(5:1)))

quintile_dat <- quintile_dat %>%
  mutate(baseline_k = remap_table$new_k[match(baseline_k, remap_table$old_k)],
         mean_k = remap_table$new_k[match(mean_k, remap_table$old_k)])

# Compare clustering assignment numbers (confusion plots) ----

all_ks <- full_join(hidim_clusts, cspline_clusts[, c("eid", "cspline_k")],
                    by = "eid")
all_ks <- full_join(all_ks, quintile_dat, by = "eid")

plotConfusionMat <- function (cmat, c1_label, c2_label) {
  nclust <- length(unique(cmat$c1))
  cmat <- cmat %>% mutate(c1 = factor(c1, levels = 1:nclust),
                          c2 = factor(c2, levels = nclust:1))
  res <- ggplot(cmat, aes(x = c1, y = c2, fill = count)) +
    geom_tile() + 
    geom_text(aes(label = count), size = 2) +
    scale_fill_gradient(low = "white", high = "#009194", guide = "none") +
    labs(x = c1_label, y = c2_label) 
  return (res)
}

ktypes <- c("hidim_k", "cspline_k", "baseline_k", "mean_k")

confusion_plots <- lapply(ktypes, function (k1type) {
  plot_list <- lapply(ktypes, function (k2type) {
    cmat <- as.data.frame(table(all_ks[, k1type],
                                all_ks[, k2type]))
    colnames(cmat) <- c("c1", "c2", "count")
    return (plotConfusionMat(cmat, c1_label = k1type, c2_label = k2type))
  })
  return (ggarrange(plotlist = plot_list, 
                    ncol = 4))
})

# Compare population observed trajectories ----

getPopnDat <- function (clust_type) {
  tmp_clust <- all_ks %>% rename(cluster = !!as.symbol(clust_type))
  full_dat <- dat %>% left_join(tmp_clust[, c("eid", "cluster")],
                                by = "eid")
  
  # Round ages to nearest 0.25 yrs
  summ_dat <- full_dat %>% 
    mutate(age_bin = plyr::round_any(age_event, 0.25, f = round)) %>%
    filter(age_bin >= 30 & age_bin <= 70) %>%
    group_by(cluster, age_bin) %>% 
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
              n = n()) %>%
    mutate(lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  # Get rolling average mean across 2 years 
  summ_dat <- summ_dat %>% 
    ungroup() %>% 
    group_by(cluster) %>% 
    arrange(age_bin, .by_group = T) %>%
    mutate(interval_width = seq_along(age_bin) - 
             findInterval(age_bin - 2, age_bin),
           mean_value_rolled = rollapply(mean_value, interval_width, mean, 
                                         fill = NA),
           cluster = factor(as.character(cluster)))
  
  summ_dat$model_type <- clust_type
  return (summ_dat)
}

for_plot <- bind_rows(getPopnDat("hidim_k"),
                      getPopnDat("cspline_k"),
                      getPopnDat("baseline_k"),
                      getPopnDat("mean_k")) 
for_plot <- for_plot %>% mutate(model_type = 
                                  factor(model_type, levels = ktypes))

popn_plot <- ggplot(for_plot,
                    aes(x = age_bin, y = mean_value, 
                        color = cluster, fill = cluster)) +
  facet_wrap(~model_type, nrow = 2, ncol = 2) +
  # Faintly plot mean and SE in background
  geom_ribbon(aes(ymin = lci_mean, 
                  ymax = uci_mean), alpha = 0.2, 
              linetype = 0) +
  geom_line(aes(y = mean_value), alpha = 0.2) +
  # Add a thicker line for rolling average 
  geom_line(aes(y = mean_value_rolled), alpha = 1) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  labs(x = "Age (years)", 
       y = paste0(PHENO, " (mean and 95 % C.I. of mean)"))

# Compare year of birth and age at first measurement in clusters ----

covars_with_k <- all_ks %>% pivot_longer(cols = -eid,
                                         names_to = "model_type",
                                         values_to = "cluster") %>%
  mutate(model_type = factor(model_type, levels = ktypes),
         cluster = factor(as.character(cluster)))
covars_with_k$baseline_age <- covars$baseline_age[match(covars_with_k$eid,
                                                        covars$eid)]
covars_with_k$year_of_birth <- general_covars$year_of_birth[match(covars_with_k$eid,
                                                                  general_covars$eid)]

covars_with_k <- covars_with_k[complete.cases(covars_with_k), ]

# Year of birth stratified by cluster identity
yob_plot <- ggplot(covars_with_k, 
                   aes(x = cluster, y = year_of_birth)) +
  facet_wrap(~model_type, nrow = 2, ncol = 2) +
  geom_violin(aes(fill = cluster), position = position_dodge(1)) +
  geom_boxplot(width = 0.1) + 
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "Cluster", y = "Birth year") 

# Baseline age stratified by cluster identity
bl_age_plot <- ggplot(covars_with_k, 
                      aes(x = cluster, y = baseline_age)) +
  facet_wrap(~model_type, nrow = 2, ncol = 2) +
  geom_violin(aes(fill = cluster), position = position_dodge(1)) +
  geom_boxplot(width = 0.1) + 
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "Cluster", y = paste0("Age ", PHENO, 
                                 " first measured (years)"))

# Print all plots ---- 

pdf(paste0(plotdir, "hidim_vs_cspline_", PHENO, "_", SEX_STRATA, ".pdf"),
    onefile = T)
print(popn_plot)
print(yob_plot)
print(bl_age_plot)
print(ggarrange(plotlist = confusion_plots, 
                nrow = 4))
dev.off()
