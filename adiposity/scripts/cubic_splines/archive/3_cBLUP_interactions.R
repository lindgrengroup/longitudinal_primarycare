# Author: Samvida S. Venkatesh
# Date: 09/11/21

library(tidyverse)
theme_set(theme_bw())

set.seed(091121)

# Read files ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1

blups <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/cubic_spline_blups_",
                 p, ".rds"))
})
names(blups) <- PHENOTYPES

SEX_STRATA <- c("F", "M", "sex_comb")

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]

# Wrangle data to add in b0 quartiles and other blup deciles ----

cut_blups <- lapply(PHENOTYPES, function (p) {
  full_dat <- lapply(SEX_STRATA, function (sx) {
    relevant_blups <- blups[[p]][[sx]]
    
    b0s <- relevant_blups[, "(Intercept)"]
    b0_quartile <- cut(b0s, 
                       breaks = c(quantile(b0s, probs = seq(0, 1, by = 1/4))),
                       labels = c(1:4),
                       include.lowest = T)
    res <- data.frame(eid = rownames(relevant_blups),
                      b0_quartile = b0_quartile)
    res <- bind_cols(res, 
                     relevant_blups[, colnames(relevant_blups) != "(Intercept)"])
    colnames(res) <- c("eid", "b0_quartile", 
                       colnames(relevant_blups)[colnames(relevant_blups) != 
                                                  "(Intercept)"])
    return (res)
  })
  names(full_dat) <- SEX_STRATA
  return (full_dat)
})
names(cut_blups) <- PHENOTYPES

# Plot raw data for selected individuals ----

## Function to get raw data for ids to plot ----

get_rawdat <- function (p, sx, b0q, spline_term, spline_d, n = 10) {
  ids_plot <- cut_blups[[p]][[sx]] %>% filter(b0_quartile == b0q)
  bs <- ids_plot[, spline_term]
  # Calculate spline term decile within this quartile of b0s
  bs_decile <- cut(bs, 
                   breaks = c(quantile(bs, 
                                       probs = seq(0, 1, by = 1/10))),
                   labels = c(1:10),
                   include.lowest = T)
  ids_plot$bs_decile <- bs_decile
  ids_plot <- ids_plot %>% filter(bs_decile == spline_d)
  
  if (n == "all") {
    ids <- ids_plot$eid
  } else {
    nsample <- min(nrow(ids_plot), n)
    ids <- sample(ids_plot$eid, nsample, replace = F)
  }
  
  res <- dat[[p]] %>% filter(eid %in% ids)
  res$b0_quartile <- b0q
  res$bs_decile <- spline_d
  res$spline_term <- spline_term
  
  return (res)
}

## Functions to create individual plots ----

plot_indivs <- function (p, sx) {
  # Get individual-level data for all quartiles, top and bottom decile
  spline_terms <- colnames(blups[[p]][[sx]])
  spline_terms <- spline_terms[spline_terms != "(Intercept)"]
  all_plots <- lapply(spline_terms, function (st) {
    res <- lapply(c(1:4), function (b0q) {
      rawdf_d1 <- get_rawdat(p, sx, b0q, st, 1, n = 10)
      rawdf_d10 <- get_rawdat(p, sx, b0q, st, 10, n = 10)
      rawdfs <- bind_rows(rawdf_d1, rawdf_d10)
      return (rawdfs)
    })
    res <- bind_rows(res)
    
    res$bs_decile <- factor(res$bs_decile)
    res$b0_quartile <- factor(res$b0_quartile)
    
    res_plot <- ggplot(res, aes(x = age_event, y = value,
                                group = eid, color = bs_decile)) +
      facet_wrap(~b0_quartile, ncol = 2) +
      geom_line() +
      geom_point() +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Age (years)", y = p, 
           title = paste0("Individuals in top or bottom deciles of ", st,  
                          " in each quartile of b0, phenotype: ",
                          p, " strata: ", sx))
    return (res_plot)
  })
  pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/cubic_spline_b0_blup_interactions_indiv_",
             p, "_", sx, ".pdf"))
  print(all_plots)
  dev.off()
}

## Function to create population-level plots ----

plot_popn <- function (p, sx) {
  # Get individual-level data for all quartiles, top and bottom decile
  spline_terms <- colnames(blups[[p]][[sx]])
  spline_terms <- spline_terms[spline_terms != "(Intercept)"]
  all_plots <- lapply(spline_terms, function (st) {
    res <- lapply(c(1:4), function (b0q) {
      rawdf_d1 <- get_rawdat(p, sx, b0q, st, 1, n = "all")
      rawdf_d10 <- get_rawdat(p, sx, b0q, st, 10, n = "all")
      rawdfs <- bind_rows(rawdf_d1, rawdf_d10)
      return (rawdfs)
    })
    res <- bind_rows(res)
    
    res$bs_decile <- factor(res$bs_decile)
    res$b0_quartile <- factor(res$b0_quartile)
    
    # Get population-level mean and S.E.M. by cutting the ages into 
    # bins of 5-year intervals for plotting 
    # Calculate mean and SE in each 5-year interval within each quartile
    res$age_bin <- cut(res$age_event, 
                       seq(20, 80, by = 5), include.lowest = T)
    summ_dat <- res %>% group_by(b0_quartile, bs_decile, age_bin) %>% 
      summarise(count = n(),
                mean_value = mean(value),
                se_value = sd(value)/sqrt(count))
    
    res_plot <- ggplot(summ_dat, aes(x = age_bin, y = mean_value, 
                                     color = bs_decile, fill = bs_decile,
                                     group = bs_decile)) +
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
                          st, " in each quartile of b0, phenotype: ",
                          p, " strata: ", sx))
    return (res_plot)
  })
  pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/cubic_splines/cubic_spline_b0_blup_interactions_popn_",
             p, "_", sx, ".pdf"))
  print(all_plots)
  dev.off()
}

## Apply plotting functions ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    plot_indivs(p, sx)
    plot_popn(p, sx)
  })
})
