# Author: Samvida S. Venkatesh
# Date: 10/03/21
# ADAPTED FROM ADIPOSITY CHANGE CONSORTIUM SOP

library(tidyverse)
theme_set(theme_bw())

# Read stratified covariate data and adjustment models ----

covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/raw_slopes_and_covars.rds")
models <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adj_slope_models.rds")

PHENOTYPES <- names(models)
STRATA <- names(models[[1]])

# Add slope residuals to covariates dataframe ----

slope_resids <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    ms <- models[[p]][[s]]
    covs <- covars[[p]][[s]]
    # Pick the residuals from the most suitable model
    res <- covs
    res$residual <- ms$m4a3$residuals
    # Flag the gainers 
    res$gainer <- res$residual > 0
    return (res)
  })
  names(res) <- STRATA
  return (res)
})
names(slope_resids) <- PHENOTYPES

# Print table of gainer stats and plot residual distributions ----

lapply(PHENOTYPES, function (p) {
  lapply(STRATA, function (s) {
    sink(paste0("log_files/adj_slopes_gainers_", p, ".txt"), append = T)
    cat(paste0("Strata: ", s, " ", "\n"))
    print(table(slope_resids[[p]][[s]]$gainer))
    cat("\n")
    sink()
  })
})

resid_plots <- lapply(PHENOTYPES, function (p) {
  df <- lapply(STRATA, function (s) {
    res <- slope_resids[[p]][[s]]
    res$sex_plot <- strsplit(s, "_")[[1]][2]
    return (res)
  })
  df <- bind_rows(df)
  p <- ggplot(df, aes(x = ancestry, y = residual)) +
    facet_wrap(~sex_plot, nrow = 3) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "adjusted slope residuals", title = p) +
    theme(legend.position = "none")
  return (p)
})

pdf("plots/adjusted_slopes/residual_distributions.pdf",
    onefile = T)
print(resid_plots)
dev.off()

saveRDS(slope_resids, "/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_adj_slopes.rds")

# Split dfs into full population, gainers, and losers ----

slope_resids <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- slope_resids[[p]][[s]]
    df_list <- list()
    df_list$all <- df
    df_list$gainers <- df[df$gainer, ]
    df_list$losers <- df[!df$gainer, ]
    return (df_list)
  })
  names(res) <- STRATA
  return (res)
})
names(slope_resids) <- PHENOTYPES
GROUPS <- c("all", "gainers", "losers")

# RINT residuals and gather final phenotypes ----

RINTed <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    res <- lapply(GROUPS, function (gp) {
      df <- slope_resids[[p]][[s]][[gp]][, c("eid", "FUyrs", "residual")]
      # Calculate FUyr quartiles
      df$FUyr_quartile <- cut(df$FUyrs, quantile(df$FUyrs), 
                              include.lowest = T, labels = 1:4)
      # Calculate Z-score within each FUyr quartile
      df <- df %>% group_by(FUyr_quartile) %>% 
        mutate(Zscore = (residual - mean(residual)) / sd(residual))
      # Apply RINT on residuals
      df$RINT_resids <- qnorm((rank(df$residual) - 0.5) / 
                                sum(!is.na(df$residual)))
      # Apply RINT on Z-scores (normalised residuals)
      df$RINT_zscores <- qnorm((rank(df$Zscore) - 0.5) / sum(!is.na(df$Zscore)))
      
      colnames(df) <- c("eid", "FUyrs", "residual", "FUyr_quartile",
                        "Zscore", "RINTed_resids", "RINTed_Zscores")
      
      write.table(df, paste0("slopes_for_analysis/", p, "_", s, "_", gp, ".txt"),
                  sep = "\t", quote = F, row.names = F)
      
      return (df)
    })
    names(res) <- GROUPS
    return (res)
  })
  names(res) <- STRATA
  return (res)
})
names(RINTed) <- PHENOTYPES

saveRDS(RINTed, "/well/lindgren/UKBIOBANK/samvida/adiposity/RINTed_final_phenotypes.rds")