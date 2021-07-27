# Author: Samvida S. Venkatesh
# Date: 02/03/21
# ADAPTED FROM ADIPOSITY CHANGE CONSORTIUM SOP

library(tidyverse)

# Read stratified covariate and median data ----

median_traits <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/covars_with_median_trait.rds")

PHENOTYPES <- names(median_traits)
STRATA <- names(median_traits[[1]])
NPCs <- 21
PCs <- paste0("PC", 1:NPCs)

median_traits <- lapply(median_traits, function (x) {
  res <- lapply(x, function (a) {
    res <- a
    res$median_age_sq <- res$median_age^2
    return (res)
  })
  return (res)
})

# Build final models ----

formula_strings <- list()
# Model covariates for BMI and WHR: 
# age and age-sq, # genetic PCs, (sex)
formula_strings$mA <- paste0("median_trait ~ median_age + median_age_sq + ", 
                             paste(PCs, collapse = " + "))
formula_strings$mAsex <- paste0("median_trait ~ median_age + median_age_sq + sex + ", 
                                paste(PCs, collapse = " + "))

# Model covariates for weight and waist circumference: 
# age and age-sq, median adult height,
# genetic PCs, (sex)
formula_strings$mB <- paste0("median_trait ~ median_age + median_age_sq + height + ", 
                             paste(PCs, collapse = " + "))
formula_strings$mBsex <- paste0("median_trait ~ median_age + median_age_sq + sex + height + ", 
                                paste(PCs, collapse = " + "))

# Run models ----

models <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- median_traits[[p]][[s]]
    
    # Model A for BMI and WHR
    if (p == "BMI" | p == "WHR") {
      # Sex-combined model
      if (grepl("sexcomb", s)) {
        run_model <- lm(formula(formula_strings$mAsex), df)
        } else { # Sex-specific models
        run_model <- lm(formula(formula_strings$mA), df)
      } 
    } else { # Model B for weight and WC
      # Sex-combined model
      if (grepl("sexcomb", s)) {
        run_model <- lm(formula(formula_strings$mBsex), df)
      } else { # Sex-specific models
        run_model <- lm(formula(formula_strings$mB), df)
      } 
    }
    return (run_model)
  })
  names(res) <- STRATA
  return (res)
})
names(models) <- PHENOTYPES

# Add model residuals to covariates dataframe ----

model_resids <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    ms <- models[[p]][[s]]
    covs <- median_traits[[p]][[s]]
    # Pick the residuals from the most suitable model
    res <- covs
    res$residual <- ms$residuals
    # Flag the gainers 
    res$gainer <- res$residual > 0
    return (res)
  })
  names(res) <- STRATA
  return (res)
})
names(model_resids) <- PHENOTYPES

# Print table of gainer stats and plot residual distributions ----

lapply(PHENOTYPES, function (p) {
  lapply(STRATA, function (s) {
    sink(paste0("log_files/adj_medians_gainers_", p, ".txt"), append = T)
    cat(paste0("Strata: ", s, " ", "\n"))
    print(table(model_resids[[p]][[s]]$gainer))
    cat("\n")
    sink()
  })
})

resid_plots <- lapply(PHENOTYPES, function (p) {
  df <- lapply(STRATA, function (s) {
    res <- model_resids[[p]][[s]]
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
    labs(x = "ancestry", y = "adjusted median residuals", title = p) +
    theme(legend.position = "none")
  return (p)
})

pdf("plots/cross_sectional/residual_distributions.pdf",
    onefile = T)
print(resid_plots)
dev.off()

saveRDS(model_resids, "/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_adj_medians.rds")

