# Author: Samvida S. Venkatesh
# Date: 18/01/22

library(tidyverse)
theme_set(theme_bw())

set.seed(180122)

# Read files ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1
SEX_STRATA <- c("F", "M", "sex_comb")

slope_models <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/without_pcs_cubic_time_splines_",
                 p, ".rds"))
})
names(slope_models) <- PHENOTYPES

blups <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/without_pcs_cubic_time_spline_blups_",
                 p, ".rds"))
})
names(blups) <- PHENOTYPES

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

COVARS_LIST <- c("baseline_age", "baseline_trait", "FUyrs", "FU_n")

# Get PC1 scores for each individual ----

pc1_scores <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    for_pca <- as.matrix(blups[[p]][[sx]])
    # Calculate PCs
    pca_res <- prcomp(for_pca, scale = F)
    pca_res <- as.data.frame(pca_res$x)
    # Get PC scores and respective eid
    to_save <- data.frame(eid = rownames(pca_res),
                     PC1_score = pca_res$PC1)
    write.table(to_save, 
                paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/PC1_scores_",
                       p, "_", sx, ".txt"),
                sep = "\t", row.names = F, quote = F)
    return (to_save)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(pc1_scores) <- PHENOTYPES

# Plot associations of PC1 scores with covariates -----

pc1_scores <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/PC1_scores_",
                      p, "_", sx, ".txt"),
               header = T, sep = "\t", stringsAsFactors = F)
    return (df)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(pc1_scores) <- PHENOTYPES

assoc_plots <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    pcs_dat <- pc1_scores[[p]][[sx]]
    covar_dat <- covars[[p]]
    plot_dat <- merge(pcs_dat, covar_dat, by = "eid")
    # Randomly select subset of ids to plot,
    # as there is otherwise too much data to plot
    nsample <- min(nrow(plot_dat), 500)
    plot_dat <- plot_dat[sample(nrow(plot_dat), nsample,
                                  replace = F), ]
    plot_dat <- plot_dat %>% 
      select(any_of(c(COVARS_LIST, "eid", "PC1_score"))) %>%
      pivot_longer(cols = all_of(COVARS_LIST), 
                   names_to = "covariate",
                   values_to = "covar_value")

    all_covars_plot <- ggplot(plot_dat, 
                          aes(x = covar_value, y = PC1_score)) +
      facet_wrap(~covariate, scales = "free_x") +
      geom_point() + 
      geom_smooth(method = "lm", formula = y~x) +
      labs(x = "Covariate value", y = "PC1 score", 
           title = paste0("PC1 score vs covariates in ",
                          sx, " phenotype: ", p))
    
    return (all_covars_plot)
  }) 
  pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/PCA/PC1_assocns_not_adj_genetic_PCs_",
             p, ".pdf"))
  print(per_sex)
  dev.off()
  return ()
})

# Plot trajectories by PC1 chunk ----

## Wrangle data ----

pc_q_membership <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    pcs_dat <- pc1_scores[[p]][[sx]] %>% 
      mutate(decile = ntile(PC1_score, 10),
             k = ifelse(decile == 10, "bottom", ifelse(decile == 1, "top",
                                                      "remove"))) %>%
      filter(k != "remove") %>% group_by(k) %>% mutate(n_cluster = n())
    covar_dat <- covars[[p]]
    plot_dat <- merge(pcs_dat, covar_dat, by = "eid")
    return (plot_dat)
  }) 
  names(per_sex) <- SEX_STRATA
  return (per_sex)
})
names(pc_q_membership) <- PHENOTYPES

##  Get random IDs within each cluster ----

get_rand_eids <- function (p, sx, n_each = 5) {
  # Return df of id and cluster number
  sampled_ids <- pc_q_membership[[p]][[sx]] %>% 
    group_by(k) %>%
    sample_n(ifelse(n_cluster < n_each, n_cluster, n_each)) %>%
    select(eid, k)
  return (sampled_ids)
}

## Function to create predicted data for set of ids ----

create_prediction_df <- function (p, sx, ids) {
  sub_dat <- subset(model_dat[[p]], 
                    model_dat[[p]]$eid %in% ids)
  # Calculate maximum time-point to predict to for each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(max_t = ceiling(max(t)))
  
  # Create new data to predict from
  new_data <- sub_dat %>% 
    select(all_of(c("eid", "baseline_age", "age_sq", "data_provider", "max_t", "sex"))) %>% 
    distinct(eid, data_provider, .keep_all = T) 
  # Timepoints to extend to
  ts <- lapply(1:nrow(new_data), FUN = function (i) { 
    seq(0, new_data$max_t[i], by = 0.25) })
  length_ts <- unlist(lapply(ts, function (x) length(x)))
  t <- unlist(ts)
  tmp <- data.frame(eid = rep(new_data$eid, length_ts),
                    data_provider = rep(new_data$data_provider, length_ts),
                    t = t)
  new_data <- merge(tmp, new_data, by = c("eid", "data_provider"))
  
  # Predict new values
  fitted_results <- 
    as.data.frame(predict(slope_models[[p]][[sx]],
                          newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # At each time-point, average across data providers for each individual
  pred_df <- pred_df %>% group_by(eid, t) %>%
    summarise(fit = mean(fit, na.rm = T)) %>%
    left_join(covars[[p]], by = "eid") %>%
    mutate(age_event = t + baseline_age) %>%
    select(eid, t, age_event, fit)
  
  return (pred_df)
}

## Function to create plots ----

plot_predictions <- function (p, sx, id_df) {
  raw_dat <- model_dat[[p]] %>% filter(eid %in% id_df$eid)
  
  plot_dat <- create_prediction_df(p, sx, id_df$eid)
  plot_dat$cluster <- id_df$k[match(plot_dat$eid,
                                    id_df$eid)]
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  
  # Order the ids by cluster number for plot
  id_levels <- id_df$eid[order(id_df$k)]
  plot_dat$eid_f <- factor(as.character(plot_dat$eid), levels = id_levels)
  raw_dat$eid_f <- factor(as.character(raw_dat$eid), levels = id_levels)
  
  res <- ggplot(plot_dat, aes(x = age_event)) +
    facet_wrap(.~eid_f, ncol = 5, scales = "free_x") +
    geom_point(data = raw_dat, aes(y = value),
               colour = "black") +
    geom_line(aes(y = fit, colour = cluster)) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Age (years)",
         y = p,
         title = paste0("Randomly selected ", sx, ", phenotype: ", p))
  return (res)
}

# Plot model predictions for hundreds of samples ----

plot_samples_modelled <- function (p, sx) {
  # Get up to 100 random IDs within each cluster
  ids_to_plot <- get_rand_eids(p, sx, 100)
  # Within each cluster, choose 5 random ids to highlight
  highlight_ids <- ids_to_plot %>% group_by(k) %>% sample_n(5)
  highlight_ids <- unique(highlight_ids$eid)
  
  # Create fitted data
  plot_dat <- create_prediction_df(p, sx, ids_to_plot$eid)
  plot_dat$cluster <- ids_to_plot$k[match(plot_dat$eid,
                                          ids_to_plot$eid)]
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  plot_dat$highlight <- plot_dat$eid %in% highlight_ids
  
  # Plot
  res_plot <- ggplot(plot_dat %>% filter(!highlight), 
                     aes(x = age_event, y = fit, group = eid,
                         colour = cluster)) +
    facet_wrap(~cluster) +
    geom_line(alpha = 0.1) +
    geom_line(data = plot_dat %>% filter(highlight),
              alpha = 1) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Age (years)",
         y = p,
         title = paste0("Sample of modelled trajectories in top/bottom PC1-quantile of phenotype: ",
                        p, " strata: ", sx))
  return (res_plot)
}

# Plot mean model prediction in each cluster ----

plot_all_modelled <- function (p, sx) {
  
  # Get list of individuals in strata
  cluster_info <- pc_q_membership[[p]][[sx]]
  clusts <- unique(cluster_info$k)
  # Get predictions in each cluster
  per_clust_pred <- lapply(clusts, function (k) {
    clust_ids <- cluster_info$eid[cluster_info$k == k]
    pred_df <- create_prediction_df(p, sx, clust_ids)
    pred_df$cluster <- k
    return (pred_df)
  })
  # Bind results across all clusters
  pred_dat <- bind_rows(per_clust_pred)
  
  # At each time-point, average across all fitted data to get mean and s.e.
  # of mean (95% CI) within each cluster
  plot_dat <- pred_dat %>% ungroup() %>%
    group_by(cluster, t) %>% summarise(mean_fit = mean(fit),
                                       sd_fit = sd(fit), 
                                       n = n()) %>%
    mutate(lci = mean_fit - 1.96*sd_fit,
           uci = mean_fit + 1.96*sd_fit,
           lci_mean = mean_fit - 1.96*(sd_fit/sqrt(n)),
           uci_mean = mean_fit + 1.96*(sd_fit/sqrt(n))) 
  # Get rolling average of these stats across 2 years 
  # because the plots are too choppy otherwise
  fn_roll <- function (x) { rollapply(x, 8, mean, fill = NA) }
  plot_dat <- plot_dat %>% ungroup() %>% 
    group_by(cluster) %>% arrange(t, .by_group = T) %>%
    mutate(across(c(mean_fit, lci, uci, lci_mean, uci_mean),
                  fn_roll))
  
  # Plot
  plot_dat$cluster <- as.factor(as.character(plot_dat$cluster))
  
  facet_plot <- ggplot(plot_dat,
                       aes(x = t, y = mean_fit, 
                           fill = cluster, 
                           colour = cluster)) +
    facet_wrap(~cluster) +
    geom_line() +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1) +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.3) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Time from baseline measurement (years)", 
         y = paste0("Mean and 95% C.I. of ", p), 
         title = paste0("Modelled trajectories in each PC1-quantile of phenotype: ",
                        p, " strata: ", sx))
  
  all_plot <- ggplot(plot_dat,
                     aes(x = t, y = mean_fit, 
                         fill = cluster, 
                         colour = cluster)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.1) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Time from baseline measurement (years)", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Modelled trajectories in each PC1-quantile of phenotype: ",
                        p, " strata: ", sx))
  
  return (list(facet_plot, all_plot))
}

# Plot population-level observed trajectories ----

plot_mean_traj_by_cluster <- function (p, sx) {
  raw_dat_with_k <- model_dat[[p]] %>%
    filter(eid %in% pc_q_membership[[p]][[sx]]$eid) 
  raw_dat_with_k$cluster <- 
    pc_q_membership[[p]][[sx]]$k[match(raw_dat_with_k$eid,
                                       pc_q_membership[[p]][[sx]]$eid)]
  # Cut ages into 1 year bins
  cut_points <- seq(20, 80, by = 1)
  cut_labels <- seq(20, 79, by = 1)
  raw_dat_with_k$age_bin <- cut(raw_dat_with_k$age_event, 
                                breaks = cut_points, 
                                include.lowest = T,
                                labels = cut_labels)
  summ_dat <- raw_dat_with_k %>% 
    group_by(cluster, age_bin) %>% 
    summarise(mean_value = mean(value),
              sd_value = sd(value), 
              n = n()) %>%
    mutate(lci = mean_value - 1.96*sd_value,
           uci = mean_value + 1.96*sd_value,
           lci_mean = mean_value - 1.96*(sd_value/sqrt(n)),
           uci_mean = mean_value + 1.96*(sd_value/sqrt(n))) 
  
  # Get rolling average of these stats across 5 years 
  # because the plots are too choppy otherwise
  fn_roll <- function (x) { rollapply(x, 5, mean, fill = NA) }
  summ_dat <- summ_dat %>% ungroup() %>% 
    group_by(cluster) %>% 
    mutate(age_bin = as.numeric(as.character(age_bin))) %>%
    arrange(age_bin, .by_group = T) %>%
    mutate(across(c(mean_value, lci, uci, lci_mean, uci_mean),
                  fn_roll))
  
  summ_dat$cluster <- as.factor(as.character(summ_dat$cluster))
  
  facet_plot <- ggplot(summ_dat, aes(x = age_bin, y = mean_value, 
                                     color = cluster, fill = cluster)) +
    facet_wrap(~cluster) +
    geom_line() +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1) +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.3) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of ", p), 
         title = paste0("Observed adiposity in each cluster of phenotype: ",
                        p, " strata: ", sx))
  
  all_plot <- ggplot(summ_dat,
                     aes(x = age_bin, y = mean_value, 
                         color = cluster, fill = cluster)) +
    geom_line() +
    geom_ribbon(aes(ymin = lci_mean, ymax = uci_mean), alpha = 0.1) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Age (years)", 
         y = paste0("Mean and 95% C.I. of mean of ", p), 
         title = paste0("Observed trajectories in top/bottom PC1-quantile of phenotype: ",
                        p, " strata: ", sx))
  
  return (list(facet_plot, all_plot))
  
}

plot_list <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(SEX_STRATA, function (sx) {
    
    # Modelled and observed trajectories of random ids
    rand_plots <- plot_predictions(p, sx, get_rand_eids(p, sx, 10))
    
    # Modelled trajectories of random ids by cluster
    cluster_sample_plots <- plot_samples_modelled(p, sx)
    
    # Modelled trajectories (all in cluster)
    modelled_all_plots <- plot_all_modelled(p, sx)
    
    # Population-level observed trajectories
    popn_plots <- plot_mean_traj_by_cluster(p, sx)
    
  
    pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/PCA/PC1_trajectories_",
               p, "_", sx, ".pdf"))
    print(rand_plots)
    print(cluster_sample_plots)
    print(modelled_all_plots)
    print(popn_plots)
    dev.off()
    
    return ()
  })
  return ()
})
