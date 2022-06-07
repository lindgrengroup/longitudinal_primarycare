# Author: Samvida S. Venkatesh
# Date: 15/11/2021

library(tidyverse)
library(lme4)
library(splines)
theme_set(theme_bw())

set.seed(151121)

# Read files ----

PHENOTYPES <- c("BMI", "Weight") 

slope_models <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_splines_",
                 p, ".rds"))
})
names(slope_models) <- PHENOTYPES

blups <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_spline_blups_",
                 p, ".rds"))
})
names(blups) <- PHENOTYPES

SEX_STRATA <- c("F", "M", "sex_comb")

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

PCs <- paste0("PC", 1:21)
COVARS <- c("baseline_age", "age_sq", "data_provider")

# Calculate PCs ----

pca_results <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SEX_STRATA, function (sx) {
    relevant_dat <- blups[[p]][[sx]] %>% select(-any_of(c("(Intercept)")))
    for_pca <- as.matrix(relevant_dat)
    # Calculate PCs
    pca_res <- prcomp(for_pca, scale = F)
    # Get PC coordinates for each individual
    pca_indiv <- pca_res$x
    return (pca_indiv)
  })
  names(res) <- SEX_STRATA
  return (res)
})
names(pca_results) <- PHENOTYPES

saveRDS(pca_results, 
        "/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_spline_blup_PCs.rds")

# Plot model predictions for individuals with high and low PC coordinates ----

# Add in covariates
model_dat <- lapply(PHENOTYPES, function (p) {
  res <- merge(dat[[p]], covars[[p]], by = "eid")
  return (res)
})
names(model_dat) <- PHENOTYPES

get_tail_PC_ids <- function (p, sx, pc, high = T, n = 25) {
  # Sample from top or bottom 10% of PC
  relevant_dat <- data.frame(eid = rownames(pca_results[[p]][[sx]]),
                             pc_value = pca_results[[p]][[sx]][, pc])
  if (high) { 
    threshold <- quantile(relevant_dat$pc_value, 0.9) 
    relevant_eids <- relevant_dat$eid[relevant_dat$pc_value > threshold]
  } else { 
    threshold <- quantile(relevant_dat$pc_value, 0.1) 
    relevant_eids <- relevant_dat$eid[relevant_dat$pc_value < threshold]
  }
  nsample <- min(length(relevant_eids), n)
  sample_eids <- sample(relevant_eids, nsample, replace = F)
  return (sample_eids)
}

get_comb_PC_ids <- function (p, sx, high_pc1 = T, high_pc2 = T, n = 9) {
  # Sample individuals with high PC1+low PC2, high PC1+high PC2, etc.
  relevant_dat <- as.data.frame(pca_results[[p]][[sx]])
  relevant_dat$eid <- rownames(pca_results[[p]][[sx]])
  
  if (high_pc1 & high_pc2) { 
    t1 <- quantile(relevant_dat$PC1, 0.9)
    t2 <- quantile(relevant_dat$PC2, 0.9)
    relevant_eids <- relevant_dat$eid[relevant_dat$PC1 > t1 & 
                                        relevant_dat$PC2 > t2]
  } else if (high_pc1 & !high_pc2) { 
    t1 <- quantile(relevant_dat$PC1, 0.9)
    t2 <- quantile(relevant_dat$PC2, 0.1)
    relevant_eids <- relevant_dat$eid[relevant_dat$PC1 > t1 &
                                        relevant_dat$PC2 < t2]
  } else if (!high_pc1 & high_pc2) {
    t1 <- quantile(relevant_dat$PC1, 0.1)
    t2 <- quantile(relevant_dat$PC2, 0.9)
    relevant_eids <- relevant_dat$eid[relevant_dat$PC1 < t1 &
                                        relevant_dat$PC2 > t2]
  } else {
    t1 <- quantile(relevant_dat$PC1, 0.1)
    t2 <- quantile(relevant_dat$PC2, 0.1)
    relevant_eids <- relevant_dat$eid[relevant_dat$PC1 < t1 &
                                        relevant_dat$PC2 < t2]
  }
  nsample <- min(length(relevant_eids), n)
  sample_eids <- sample(relevant_eids, nsample, replace = F)
  return (sample_eids)
}

## Function to create predicted data for set of ids ----

create_prediction_df <- function (p, sx, ids) {
  sub_dat <- subset(model_dat[[p]], 
                    model_dat[[p]]$eid %in% ids)
  # Calculate maximum time-point to predict to for each eid
  sub_dat <- sub_dat %>% group_by(eid) %>% 
    mutate(min_age = min(age_event),
           max_age = max(age_event))
  
  # Create new data to predict from
  new_data <- sub_dat %>% select(all_of(c("eid", COVARS, PCs,
                                          "min_age", "max_age", "sex"))) %>% 
    distinct(eid, data_provider, .keep_all = T) 
  # Timepoints to extend to 
  t <- unlist(lapply(1:nrow(new_data), FUN = function (i) { 
    seq(new_data$min_age[i], new_data$max_age[i], length.out = 30) }))
  tmp <- data.frame(eid = rep(new_data$eid, each = 30),
                    data_provider = rep(new_data$data_provider, each = 30),
                    age_event = t)
  new_data <- merge(tmp, new_data, by = c("eid", "data_provider"))
  
  # Predict new values
  fitted_results <- 
    as.data.frame(predict(slope_models[[p]][[sx]],
                          newdata = new_data))
  colnames(fitted_results) <- "fit"
  pred_df <- bind_cols(new_data, fitted_results)
  
  # At each time-point, average across data providers for each individual
  pred_df <- pred_df %>% group_by(eid, age_event) %>%
    summarise(fit = mean(fit))
  
  return (pred_df)
}

## Function to create plots ----

sex_col_palette <- c("#F8766D", "#00BFC4", "#C77CFF")
names(sex_col_palette) <- c("F", "M", "sex_comb")

plot_predictions <- function (p, sx, ids, numcol = 5) {
  raw_dat <- model_dat[[p]] %>% filter(eid %in% ids)
  
  plot_dat <- create_prediction_df(p, sx, ids) %>%
    mutate(model_strata = sx)
  
  res <- ggplot(plot_dat, aes(x = age_event)) +
    facet_wrap(~eid, ncol = numcol) +
    geom_point(data = raw_dat, aes(y = value),
               colour = "black") +
    geom_line(aes(y = fit, colour = model_strata)) +
    scale_color_manual(values = sex_col_palette) +
    labs(x = "Age (years)",
         y = p)
  return (res)
}

## Apply plotting functions ----

plot_list <- lapply(PHENOTYPES, function (p) {
  per_sex <- lapply(c("F", "M", "sex_comb"), function (sx) {
    
    # Low or high PCs
    pc_names <- colnames(pca_results[[p]][[sx]])
    
    low_pc_plots <- lapply(pc_names, function (pc_name) {
      res_plot <- plot_predictions(p, sx, 
                                   get_tail_PC_ids(p, sx, pc_name, 
                                                   high = F, 25),
                                   numcol = 5)
      res_plot <- res_plot +
        labs(title = paste0(sx, " in bottom 10% of ", pc_name,
                            " phenotype: ", p))
      return (res_plot)
    })
    
    high_pc_plots <- lapply(pc_names, function (pc_name) {
      res_plot <- plot_predictions(p, sx, 
                                   get_tail_PC_ids(p, sx, pc_name, 
                                                   high = T, 25),
                                   numcol = 5)
      res_plot <- res_plot +
        labs(title = paste0(sx, " in top 10% of ", pc_name,
                            " phenotype: ", p))
      return (res_plot)
    })
    
    # Combination of PCs 
    top_top <- plot_predictions(p, sx, 
                                get_comb_PC_ids(p, sx, 
                                                high_pc1 = T, high_pc2 = T, 9),
                                numcol = 3) +
      labs(title = paste0(sx, " in top 10% of PC1 and PC2, phenotype: ", p))
    
    top_bot <- plot_predictions(p, sx, 
                                get_comb_PC_ids(p, sx, 
                                                high_pc1 = T, high_pc2 = F, 9),
                                numcol = 3) +
      labs(title = 
             paste0(sx, 
                    " in top 10% of PC1 and bottom 10% of PC2, phenotype: ", 
                    p))
    
    bot_top <- plot_predictions(p, sx, 
                                get_comb_PC_ids(p, sx, 
                                                high_pc1 = F, high_pc2 = T, 9),
                                numcol = 3) +
      labs(title = 
             paste0(sx, 
                    " in bottom 10% of PC1 and top 10% of PC2, phenotype: ", 
                    p))
    
    bot_bot <- plot_predictions(p, sx, 
                                get_comb_PC_ids(p, sx, 
                                                high_pc1 = F, high_pc2 = F, 9),
                                numcol = 3) +
      labs(title = 
             paste0(sx, 
                    " in bottom 10% of PC1 and bottom 10% of PC2, phenotype: ", 
                    p))
    
    pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/cubic_spline_PCA_tails_",
               p, "_", sx, ".pdf"))
    print(low_pc_plots)
    print(high_pc_plots)
    print(top_top)
    print(top_bot)
    print(bot_top)
    print(bot_bot)
    dev.off()
    
    return ()
  })
  return ()
})

# Plot raw data for indivudals with high and low PC coordinates ----

# Wrangle data to add in b0 quartiles and PC values ----

cut_blups <- lapply(PHENOTYPES, function (p) {
  full_dat <- lapply(SEX_STRATA, function (sx) {
    # Get b0 quartile
    relevant_blups <- blups[[p]][[sx]]
    b0s <- relevant_blups[, "(Intercept)"]
    b0_quartile <- cut(b0s, 
                       breaks = c(quantile(b0s, probs = seq(0, 1, by = 1/4))),
                       labels = c(1:4),
                       include.lowest = T)
    res <- data.frame(eid = rownames(relevant_blups),
                      b0_quartile = b0_quartile)
    
    # Add in PC1 and PC2 values
    relevant_PCs <- as.data.frame(pca_results[[p]][[sx]])
    relevant_PCs$eid <- rownames(relevant_PCs)
    
    res <- merge(res, relevant_PCs, by = "eid")
    return (res)
  })
  names(full_dat) <- SEX_STRATA
  return (full_dat)
})
names(cut_blups) <- PHENOTYPES

## Function to get raw data for ids to plot ----

get_rawdat <- function (p, sx, b0q, pc, pc_dec) {
  ids_plot <- cut_blups[[p]][[sx]] %>% filter(b0_quartile == b0q)
  pc_vals <- ids_plot[, pc]
  # Calculate PC decile within this quartile of b0s
  pc_decile <- cut(pc_vals, 
                   breaks = c(quantile(pc_vals, 
                                       probs = seq(0, 1, by = 1/10))),
                   labels = c(1:10),
                   include.lowest = T)
  ids_plot$pc_decile <- pc_decile
  ids_plot <- ids_plot %>% filter(pc_decile == pc_dec)
  
  res <- dat[[p]] %>% filter(eid %in% ids_plot$eid)
  res$b0_quartile <- b0q
  res$pc_decile <- pc_dec
  res$pc_term <- pc
  return (res)
}

## Function to plot individual-level plots ----

plot_indiv <- function (p, sx) {
  # Get individual-level data for all quartiles, top and bottom decile
  pcs <- colnames(pca_results[[p]][[sx]])
  all_plots <- lapply(pcs, function (pc_term) {
    res <- lapply(c(1:4), function (b0q) {
      rawdf_d1 <- get_rawdat(p, sx, b0q, pc_term, 1)
      rawdf_d10 <- get_rawdat(p, sx, b0q, pc_term, 10)
      rawdfs <- bind_rows(rawdf_d1, rawdf_d10)
      return (rawdfs)
    })
    res <- bind_rows(res)
    
    res$pc_decile <- factor(res$pc_decile)
    res$b0_quartile <- factor(res$b0_quartile)
    
    # Within each b0 quartile and PC decile, choose 5 IDs to highlight
    hids <- res %>% group_by(pc_decile, b0_quartile) %>% sample_n(5)
    res$highlight <- res$eid %in% hids$eid
    
    alpha_scale <- c(0.1, 1)
    names(alpha_scale) <- c(F, T)
    
    res_plot <- ggplot(data = subset(res, !res$highlight), 
                       aes(x = age_event, y = value, 
                           color = pc_decile, fill = pc_decile,
                           group = eid, alpha = highlight)) +
      facet_wrap(~b0_quartile, ncol = 2) +
      geom_point() +
      geom_path() +
      geom_point(data = subset(res, res$highlight)) +
      geom_path(data = subset(res, res$highlight)) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      scale_alpha_manual(values = alpha_scale) +
      labs(x = "Age (years)", 
           y = p, 
           title = paste0("Raw values in top or bottom deciles of ",
                          pc_term, " in each quartile of b0, phenotype: ",
                          p, " strata: ", sx))
    png(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/cubic_spline_b0_",
               pc_term, "_interactions_", p, "_", sx, ".png"))
    print(res_plot)
    dev.off()
    return (res_plot)
  })
}

## Function to plot population-level plots ----

plot_popn <- function (p, sx) {
  # Get individual-level data for all quartiles, top and bottom decile
  pcs <- colnames(pca_results[[p]][[sx]])
  all_plots <- lapply(pcs, function (pc_term) {
    res <- lapply(c(1:4), function (b0q) {
      rawdf_d1 <- get_rawdat(p, sx, b0q, pc_term, 1)
      rawdf_d10 <- get_rawdat(p, sx, b0q, pc_term, 10)
      rawdfs <- bind_rows(rawdf_d1, rawdf_d10)
      return (rawdfs)
    })
    res <- bind_rows(res)
    
    res$pc_decile <- factor(res$pc_decile)
    res$b0_quartile <- factor(res$b0_quartile)
    
    # Get population-level mean and S.E.M. by cutting the ages into 
    # bins of 5-year intervals for plotting 
    # Calculate mean and SE in each 5-year interval within each quartile
    res$age_bin <- cut(res$age_event, 
                       seq(20, 80, by = 5), include.lowest = T)
    summ_dat <- res %>% group_by(b0_quartile, pc_decile, age_bin) %>% 
      summarise(count = n(),
                mean_value = mean(value),
                se_value = sd(value)/sqrt(count))
    
    res_plot <- ggplot(summ_dat, aes(x = age_bin, y = mean_value, 
                                     color = pc_decile, fill = pc_decile,
                                     group = pc_decile)) +
      facet_wrap(~b0_quartile, ncol = 2) +
      geom_point() +
      geom_path() +
      geom_ribbon(aes(ymin = mean_value - 1.96*se_value, 
                      ymax = mean_value + 1.96*se_value),
                  alpha = 0.2) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Age bin (years)", 
           y = paste0("Mean and 95% C.I. of mean of ", p), 
           title = paste0("Mean adiposity in top or bottom deciles of ",
                          pc_term, " in each quartile of b0, phenotype: ",
                          p, " strata: ", sx))
    return (res_plot)
  })
  pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/cubic_spline_b0_PC_interactions_popn_",
             p, "_", sx, ".pdf"))
  print(all_plots)
  dev.off()
}

## Apply plotting functions ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    plot_indiv(p, sx)
    plot_popn(p, sx)
  })
})

