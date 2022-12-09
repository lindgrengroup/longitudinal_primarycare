# Author: Samvida S. Venkatesh
# Date: 22/03/21

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

# Read models, covariates, and adiposity data ----

covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/raw_slopes_and_covars.rds")
models <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adj_slope_models.rds")
adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/stratified_adiposity.rds")

PHENOTYPES <- names(models)
STRATA <- names(models[[1]])

# Re-create full adiposity (non-stratified) df for plots
long_adiposity <- lapply(adiposity, function (x) bind_rows(x) )
names(long_adiposity) <- PHENOTYPES

# Add slope residuals to covariates dataframe ----

slope_resids <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    ms <- models[[p]][[s]]
    covs <- covars[[p]][[s]]
    res <- cbind.data.frame(covs, 
                              resid1 = ms$m4$residuals,
                              resid2 = ms$m4a1$residuals,
                              resid3 = ms$m4a2$residuals,
                              resid4 = ms$m4a3$residuals)
    # Long format for plotting
    res <- res %>% pivot_longer(starts_with("resid"), 
                                names_to = "model",
                                names_prefix = "resid",
                                values_to = "residual")
    res$gainer <- res$residual > 0
    return (res)
  })
  names(res) <- STRATA
  return (res)
})
names(slope_resids) <- PHENOTYPES

# Plot mean (binned) trajectories by adj slope quartile ----

as_quartile_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    rs <- slope_resids[[p]][[s]]
    rs <- rs %>% group_by(model) %>% 
      mutate(q = cut(residual, quantile(residual), include.lowest = T,
                     labels = paste0("q", 1:4)))
    
    df <- long_adiposity[[p]]
    df <- inner_join(df, rs[, c("eid", "model", "q")], by = "eid")
    
    # Calculate mean and SE in each 5-year interval within each quartile
    age_bin_cuts <- seq(20, 80, by = 5)
    df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
    plot_df <- df %>% group_by(model, q, age_bin) %>% 
      summarise(count = n(),
                mean_value = mean(value),
                se_value = sd(value)/sqrt(count))
    
    # Plot 
    p <- ggplot(plot_df, aes(x = age_bin, y = mean_value,
                             group = q, color = q, fill = q)) +
      facet_wrap(~model, nrow = 2) + 
      geom_point() +
      geom_path() +
      geom_ribbon(aes(ymin = mean_value - se_value, 
                      ymax = mean_value + se_value),
                  alpha = 0.2) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
           title = s)
    return (p)
  })
  names(res) <- STRATA
  pdf(paste0("plots/adjusted_slopes/model_comparison_quartile_status_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(as_quartile_trajectories) <- PHENOTYPES

# Plot mean (binned) trajectories by gainer status ----

gainer_status_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    rs <- slope_resids[[p]][[s]]
    df <- long_adiposity[[p]]
    
    df <- inner_join(df, rs[, c("eid", "model", "gainer")], by = "eid")
    
    # Calculate mean and SE in each 5-year interval within gainers/non-gainers
    age_bin_cuts <- seq(20, 80, by = 5)
    df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
    plot_df <- df %>% group_by(model, gainer, age_bin) %>% 
      summarise(count = n(),
                mean_value = mean(value),
                se_value = sd(value)/sqrt(count))
    
    # Plot 
    p <- ggplot(plot_df, aes(x = age_bin, y = mean_value,
                       group = gainer, color = gainer, fill = gainer)) +
      facet_wrap(~model, nrow = 2) +
      geom_point() +
      geom_path() +
      geom_ribbon(aes(ymin = mean_value - se_value, 
                      ymax = mean_value + se_value),
                  alpha = 0.2) +
      scale_color_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C"), 
                         guide = F) +
      scale_fill_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C"), 
                        guide = F) +
      labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
           title = s)
    return (p)
  })
  names(res) <- STRATA
  pdf(paste0("plots/adjusted_slopes/model_comparison_gainer_status_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(gainer_status_trajectories) <- PHENOTYPES

