# Author: Samvida S. Venkatesh
# Date: 24/05/2022

library(argparse)
library(splines)
library(tidyverse)
library(zoo)
library(pheatmap)
library(ggpubr)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

# Get arguments ----

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

# High-dimensional spline modelling results
model_dat <- readRDS(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/results/fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))
B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

# Cluster assignments
clust_res <- read.table(paste0(resdir, "assigned_clusters_", PHENO, "_", SEX_STRATA, 
                               ".txt"),
                        sep = "\t", header = T, stringsAsFactors = F)

# Data (adjusted and baselined) and raw weight data
raw_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/data/dat_to_model.rds")[[PHENO]][[SEX_STRATA]]
orig_popn_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[[PHENO]]

# Covariates
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[[PHENO]]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)

# Wrangle data ----

covars$eid <- as.character(covars$eid)
general_covars$eid <- as.character(general_covars$eid)

covars <- merge(covars, general_covars, by = "eid")

# Add in covariates to raw GP data
raw_dat$eid <- as.character(raw_dat$eid)
add_covs <- covars %>% select(any_of(c("eid", "baseline_age", 
                                       "sex", "year_of_birth", "smoking_status",
                                       "age_at_death")))
full_dat <- merge(raw_dat, add_covs, by = "eid")

if (SEX_STRATA != "sex_comb") full_dat <- full_dat %>% filter(sex == SEX_STRATA)

# Count number of individuals in each cluster to make sampling quicker later
clust_res <- clust_res %>% 
  group_by(clust) %>% 
  mutate(n_cluster = n())

# Plot raw data and model predictions for random individuals ----

## Get random IDs within each cluster
getRandIDs <- function (n_each = 5) {
  # Return df of id and cluster number
  sampled_ids <- clust_res %>% 
    group_by(clust) %>%
    sample_n(ifelse(n_cluster < n_each, n_cluster, n_each)) %>%
    mutate(eid = as.character(eid)) %>%
    select(eid, clust)
  return (sampled_ids)
}

## Create predicted data for set of ids
createPredDat <- function (id_list) {
  res_list <- lapply(id_list, function (id) {
    pred_value <- B %*% spline_posteriors[[id]]$mu
    sd_pred <- sqrt(diag(B %*% spline_posteriors[[id]]$Sig %*% t(B)) * model_resid_var)
    
    res <- data.frame(eid = id,
                      t_diff = 1:length(pred_value),
                      fit_mean = pred_value,
                      fit_sd = sd_pred)
  })
  res_list <- bind_rows(res_list)
  return (res_list)
}

## Create plots
plotPreds <- function (id_df) {
  sample_raw <- raw_dat %>% filter(eid %in% id_df$eid)
  sample_pred <- createPredDat(id_df$eid)
  
  sample_pred$clust <- id_df$clust[match(sample_pred$eid,
                                    id_df$eid)]
  sample_pred$clust <- factor(as.character(sample_pred$clust))
  
  # Order the ids by cluster number for plot
  id_levels <- id_df$eid[order(id_df$clust)]
  sample_pred$eid_f <- factor(as.character(sample_pred$eid), 
                              levels = id_levels)
  sample_raw$eid_f <- factor(as.character(sample_raw$eid), 
                          levels = id_levels)
  
  res <- ggplot(sample_pred, aes(x = t_diff)) +
    facet_wrap(~eid_f, nrow = 5, ncol = 5, scales = "free_y") +
    geom_point(data = sample_raw,
               aes(x = t_diff, y = value_fulladj)) +
    geom_line(aes(y = fit_mean,
                  colour = clust)) +
    geom_ribbon(aes(ymin = fit_mean - 1.96*fit_sd,
                    ymax = fit_mean + 1.96*fit_sd,
                    fill = clust, colour = clust), alpha = 0.1) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    scale_colour_brewer(palette = "Set1", guide = "none") +
    labs(x = "Days from first measurement", y = "Confounder-adj value")
  
  return (res)
}

# Plot cluster-level population observed trajectories ----

plotPopnTrajCluster <- function () {
  pop_dat_with_k <- orig_popn_dat
  pop_dat_with_k$cluster <- clust_res$clust[match(pop_dat_with_k$eid,
                                                  clust_res$eid)]
  # Round ages to nearest 0.25 yrs
  summ_dat <- pop_dat_with_k %>% 
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
  
  all_plot <- ggplot(summ_dat,
                     aes(x = age_bin, y = mean_value, 
                         color = cluster, fill = cluster)) +
    # Faintly plot mean and SE in background
    geom_ribbon(aes(ymin = lci_mean, 
                    ymax = uci_mean), alpha = 0.2, 
                linetype = 0) +
    geom_line(aes(y = mean_value), alpha = 0.2) +
    # Add a thicker line for rolling average 
    geom_line(aes(y = mean_value_rolled), alpha = 1) +
    scale_color_brewer(palette = "Set1", guide = F) +
    scale_fill_brewer(palette = "Set1", guide = F) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 8))
  
  return (all_plot)
}

# Plot cluster identity vs covariates ----

sex_col_palette <- c("#F8766D", "#00BFC4", "#C77CFF")
names(sex_col_palette) <- c("F", "M", "sex_comb")

plotClusterCovs <- function () {
  # Get covariates data
  covars_with_k <- covars %>%
    filter(eid %in% clust_res$eid) 
  covars_with_k$cluster <- 
    clust_res$clust[match(covars_with_k$eid,
                            clust_res$eid)]
  covars_with_k$cluster <- as.factor(as.character(covars_with_k$cluster))
  
  # Number of men/women (individuals) in each cluster
  summ_sex <- covars_with_k %>% count(cluster, sex)
  sex_plot <- ggplot(summ_sex, 
                     aes(x = cluster, y = n, fill = sex)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = sex_col_palette) +
    labs(x = "Cluster", y = "Number of individuals") 
  
  # Year of birth stratified by cluster identity
  yob_plot <- ggplot(covars_with_k, 
                        aes(x = cluster, y = year_of_birth)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Cluster", y = "Birth year") +
    theme(legend.position = "none")
  
  # Age at death stratified by cluster identity
  death_plot <- ggplot(covars_with_k, 
                     aes(x = cluster, y = age_at_death)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Cluster", y = "Age at death (years)") +
    theme(legend.position = "none")
  
  # Baseline age stratified by cluster identity
  bl_age_plot <- ggplot(covars_with_k, 
                        aes(x = cluster, y = baseline_age)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Cluster", y = paste0("Age ", PHENO, 
                                            " first measured (years)")) +
    theme(legend.position = "none")
  
  # Baseline trait value
  bl_trait_plot <- ggplot(covars_with_k, 
                          aes(x = cluster, y = baseline_trait)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Cluster", y = paste0("Baseline ", PHENO)) +
    theme(legend.position = "none")
  
  # Number of follow-up measures
  max_y <- min(50, max(covars_with_k$FU_n))
  fu_n_plot <- ggplot(covars_with_k, 
                      aes(x = cluster, y = FU_n)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    scale_y_continuous(limits = c(0, max_y)) +
    labs(x = "Cluster", y = "# follow-up measures") +
    theme(legend.position = "none")
  
  # Number of follow-up years
  fu_yrs_plot <- ggplot(covars_with_k, 
                        aes(x = cluster, y = FUyrs)) +
    geom_violin(aes(fill = cluster), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Set1") + 
    labs(x = "Cluster", y = "# years follow-up") +
    theme(legend.position = "none")
  
  # Smoking status in each cluster
  summ_smoking <- covars_with_k %>% 
    group_by(cluster, smoking_status) %>%
    summarise(n = n()) %>% 
    mutate(freq = n/sum(n))
  smoker_plot <- ggplot(summ_smoking, 
                     aes(x = cluster, y = freq, fill = smoking_status)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Cluster", y = "Proportion of individuals") 
  
  
  arranged_plots <- ggarrange(plotlist = list(sex_plot, yob_plot, death_plot,
                                              bl_age_plot, bl_trait_plot, 
                                              fu_n_plot, fu_yrs_plot,
                                              smoker_plot),
                              nrow = 2, ncol = 2)
  
  return (arranged_plots)
}

# Apply all plotting functions ----

# Modelled and observed trajectories of random ids
rand_plots <- plotPreds(getRandIDs(5))

# Population-level observed trajectories
popn_plots <- plotPopnTrajCluster()

# Cluster membership vs covariates 
cov_plots <- plotClusterCovs()

pdf(paste0(plotdir, "clustering_results_", 
           PHENO, "_", SEX_STRATA, ".pdf"), onefile = T)
print(rand_plots)
print(popn_plots)  
print(cov_plots)
dev.off()

# Plot heatmaps of spline coefficients in the clusters ----

# RINT coefficient for plotting
rint_x <- function (x) {
  return (qnorm((rank(x) - 0.5) / sum(!is.na(x))))
}

# Row annotations (cluster membership)
k_info <- clust_res
# Order b-coefs by cluster membership
k_info <- k_info[order(k_info$clust), ]
row_annots <- data.frame(cluster = as.character(k_info$clust))
rownames(row_annots) <- as.character(k_info$eid)

mn_mat <- lapply(spline_posteriors, function (spobj) {
  return (as.data.frame(t(spobj$mu)))
})
mn_mat <- bind_rows(mn_mat)
rownames(mn_mat) <- names(spline_posteriors)
# Arrange by id
mn_mat <- as.matrix(mn_mat[rownames(row_annots), ])

# For plotting, scale each column (RINT) to see differences better
mn_mat <- apply(mn_mat, 2, FUN = rint_x)

png(paste0(plotdir, "heatmaps_", PHENO, "_", SEX_STRATA, ".png"),
    width = 21, height = 29.7, units = "cm", res = 300)
pheatmap(mn_mat, 
         cluster_rows = F, cluster_cols = F,
         scale = "none",   
         annotation_row = row_annots, 
         treeheight_col = 0, treeheight_row = 50,
         show_rownames = F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         legend = T, fontsize = 10,
         main = paste0("B-coefs in", PHENO, "_", SEX_STRATA))
dev.off()

