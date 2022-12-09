# Author: Samvida S. Venkatesh
# Date: 09/11/21

library(tidyverse)
theme_set(theme_bw())

set.seed(091121)

# Read files ----

PHENOTYPES <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/pheno_names.txt")$V1

blups <- lapply(PHENOTYPES, function (p) {
  readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/results/lmm_blups_",
                 p, ".rds"))
})
names(blups) <- PHENOTYPES

SEX_STRATA <- c("F", "M", "sex_comb")

dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]

# Wrangle data to add in b0 quartiles and b1 deciles ----

cut_blups <- lapply(PHENOTYPES, function (p) {
  full_dat <- lapply(SEX_STRATA, function (sx) {
    relevant_blups <- blups[[p]][[sx]]
   
    b0s <- relevant_blups[, "(Intercept)"]
    b0_quartile <- cut(b0s, 
                       breaks = c(quantile(b0s, probs = seq(0, 1, by = 1/4))),
                       labels = c(1:4),
                       include.lowest = T)
    res <- data.frame(eid = rownames(relevant_blups),
                      b0_quartile = b0_quartile,
                      b1 = relevant_blups[, "t"])
    return (res)
  })
  names(full_dat) <- SEX_STRATA
  return (full_dat)
})
names(cut_blups) <- PHENOTYPES

# Plot raw data for selected individuals ----

## Function to get raw data for ids to plot ----

get_rawdat <- function (p, sx, b0q, b1d, n = 10) {
  ids_plot <- cut_blups[[p]][[sx]] %>% filter(b0_quartile == b0q)
  
  # Calculate slope decile within this quartile of b0s
  b1_decile <- cut(ids_plot$b1, 
                   breaks = c(quantile(ids_plot$b1, 
                                       probs = seq(0, 1, by = 1/10))),
                   labels = c(1:10),
                   include.lowest = T)
  ids_plot$b1_decile <- b1_decile
  ids_plot <- ids_plot %>% filter(b1_decile == b1d)
  
  if (n == "all") {
    ids <- ids_plot$eid
  } else {
    nsample <- min(nrow(ids_plot), n)
    ids <- sample(ids_plot$eid, nsample, replace = F)
  }
  
  res <- dat[[p]] %>% filter(eid %in% ids)
  res$b0_quartile <- b0q
  res$b1_decile <- b1d
    
  return (res)
}

## Functions to create plots ----

plot_indivs <- function (p, sx) {
  # Get individual-level data for all quartiles, top and bottom decile
  full_dat <- lapply(1:4, function (b0q) {
    res <- lapply(c(1, 10), function (b1d) {
      rawdf <- get_rawdat(p, sx, b0q, b1d, n = 10)
    })
    res <- bind_rows(res)
    return (res)
  })
  full_dat <- bind_rows(full_dat)
  
  full_dat$b1_decile <- factor(full_dat$b1_decile)
  full_dat$b0_quartile <- factor(full_dat$b0_quartile)
  
  res_plot <- ggplot(full_dat, aes(x = age_event, y = value,
                                   group = eid, color = b1_decile)) +
    facet_wrap(~b0_quartile, ncol = 2) +
    geom_line() +
    geom_point() +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Age (years)", y = p, 
         title = paste0("Individuals in top or bottom deciles of b1 in each quartile of b0, phenotype: ",
                        p, " strata: ", sx))
  return (res_plot)
}

plot_popn <- function (p, sx) {
  # Get individual-level data for all quartiles, top and bottom decile
  full_dat <- lapply(1:4, function (b0q) {
    res <- lapply(c(1, 10), function (b1d) {
      rawdf <- get_rawdat(p, sx, b0q, b1d, n = "all")
    })
    res <- bind_rows(res)
    return (res)
  })
  full_dat <- bind_rows(full_dat)
  full_dat$b1_decile <- factor(full_dat$b1_decile)
  full_dat$b0_quartile <- factor(full_dat$b0_quartile)
  # Get population-level mean and S.E.M. by cutting the ages into 
  # bins of 5-year intervals for plotting 
  # Calculate mean and SE in each 5-year interval within each quartile
  full_dat$age_bin <- cut(full_dat$age_event, 
                          seq(20, 80, by = 5), include.lowest = T)
  summ_dat <- full_dat %>% group_by(b0_quartile, b1_decile, age_bin) %>% 
    summarise(count = n(),
              mean_value = mean(value),
              se_value = sd(value)/sqrt(count))
  
  res_plot <- ggplot(summ_dat, aes(x = age_bin, y = mean_value, 
                                   color = b1_decile, fill = b1_decile,
                                   group = b1_decile)) +
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
         title = paste0("Mean adiposity in top or bottom deciles of b1 in each quartile of b0, phenotype: ",
                        p, " strata: ", sx))
  return (res_plot)
}

## Apply plotting functions ----

indiv_plots <- lapply(PHENOTYPES, function (p) {
  per_strata <- lapply(SEX_STRATA, function (sx) {
    res <- plot_indivs(p, sx)
    return (res)
  })
  return (per_strata)
})
pdf("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/lmm_b0_b1_interaction_indiv.pdf")
print(indiv_plots)
dev.off()

popn_plots <- lapply(PHENOTYPES, function (p) {
  per_strata <- lapply(SEX_STRATA, function (sx) {
    res <- plot_popn(p, sx)
    return (res)
  })
  return (per_strata)
})
pdf("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/plots/lmm_b0_b1_interaction_popn.pdf")
print(popn_plots)
dev.off()

