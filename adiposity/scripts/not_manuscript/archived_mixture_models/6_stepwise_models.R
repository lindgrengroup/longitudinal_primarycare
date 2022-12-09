# Author: Samvida S. Venkatesh
# Date: 21/05/2021

library(tidyverse)
theme_set(theme_bw())
library(flexmix)

set.seed(210521)

# Read files ----

args <- commandArgs(trailingOnly = T)

PHENO <- args[1]
STRATUM <- args[2]
STRATUM_ANCESTRY <- sub("_.*", "", STRATUM)
STRATUM_SEX <- sub(".*_", "", STRATUM)

dat <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/adj_traits_and_covars/", 
                    PHENO, "_", STRATUM, ".rds"))

# Flexmix step through k = 1:K (10) ----

# Linear model (as splines don't converge for most non-white strata)
mod_LIN <- FLXMRlmm(adj_value ~ 0 + age_event, 
                    random = ~ 1 + age_event,
                    varFix = c(Random = F, Residual = F))

# Flexmix to step through 1:K clusters
K <- 10

runFlex <- function (mod) {
  out <- tryCatch(
    expr = {
      # Run model
      res <- suppressWarnings(stepFlexmix(.~.|eid, k = 1:K, 
                                          model = mod, 
                                          data = dat, 
                                          control = list(iter.max = 1000, minprior = 0),
                                          nrep = 5))
      message("Running step-wise increasing models...")
      res
    }, 
    error = function (cond) {
      message("Models fail to run, original error message:")
      message(cond)
      # Return value in case of error
      return (NA)
    }
  )
  return (out)
}

flex_step_LIN <- runFlex(mod_LIN)

# Plot BIC for models of increasing K ----

plot_BIC_dat <- data.frame(k = 1:K,
                           BIC = BIC(flex_step_LIN))

# Report BIC plot 
pdf(paste0("plots/mixture_models/stepwise_BIC_", PHENO, "_", STRATUM, ".pdf"),
    onefile = T)
ggplot(plot_BIC_dat, aes(x = k, y = BIC)) +
  geom_point(shape = 1, size = 3) +
  geom_line() +
  scale_x_continuous(breaks = c(1:10), labels = c(1:10)) +
  labs(title = paste(PHENO, STRATUM))
dev.off()

# Save best model according to BIC ----
best_flex <- getModel(flex_step_LIN, "BIC")

saveRDS(best_flex, 
        paste0("results/mixture_models/best_BIC_model_", PHENO, "_",
               STRATUM, ".rds"))

# Examine best model ----

PHENO <- c("BMI", "WHR")
STRATA <- c("asian_F", "asian_M", "asian_sexcomb",
            "black_F", "black_M", "black_sexcomb",
            "chinese_F", "chinese_M", "chinese_sexcomb",
            "mixed_F", "mixed_M", "mixed_sexcomb",
            "other_F", "other_M", "other_sexcomb")

raw_dat <- lapply(PHENO, function (p) {
  res_list <- lapply(STRATA, function (s) {
    fname <- paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/adj_traits_and_covars/", 
                    p, "_", s, ".rds")
    if (file.exists(fname)) {
      res <- readRDS(fname)
    } else { res  <- NA }
    return (res)
  })
  names(res_list) <- STRATA
  return (res_list)
})
names(raw_dat) <- PHENO

best_flex <- lapply(PHENO, function (p) {
  res_list <- lapply(STRATA, function (s) {
    fname <- paste0("results/mixture_models/best_BIC_model_", 
                    p, "_", s, ".rds")
    if (file.exists(fname)) {
      res <- readRDS(fname)
    } else { res  <- NA }
    return (res)
  })
  names(res_list) <- STRATA
  return (res_list)
})
names(best_flex) <- PHENO

# Fitted values at group level ----

mod_preds <- function (dat, mod) {
  pred_fits <- predict(mod, newdata = dat, 
                       se.fit = T)
  pred_fits <- bind_rows(pred_fits)
  pred_fits$age_event <- dat$age_event
  res <- pivot_longer(pred_fits, cols = starts_with("Comp"),
                      names_to = "cluster", values_to = "mean")
  group_sd <- sqrt(parameters(mod)["sigma2.Residual", ])
  res <- res %>% mutate(lci = mean - group_sd[cluster],
                       uci = mean + group_sd[cluster])
  return (res)
}

pred_plots <- lapply(PHENO, function (p) {
  res_plots <- lapply(STRATA, function (s) {
    if (!is.na(best_flex[[p]][[s]])) {
      rdat <- raw_dat[[p]][[s]]
      pred_dat <- data.frame(age_event = seq(min(rdat$age_event), 
                                             max(rdat$age_event), 
                                             length.out = 100),
                             eid = 1)
      pred_dat <- mod_preds(pred_dat, best_flex[[p]][[s]])
      
      p <- ggplot(pred_dat, aes(x = age_event, y = mean, 
                                colour = cluster, fill = cluster)) +
        geom_line() +
        geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
        scale_colour_brewer(palette = "Set1") + 
        scale_fill_brewer(palette = "Set1") + 
        labs(x = "Age", y = "Predicted adj. adiposity", title = s)
    } else { p <- NA }
    return (p)
  })
  pdf(paste0("plots/mixture_models/best_model_fits_", p, ".pdf"), 
      onefile = T)
  print(res_plots)
  dev.off()
  return (res_plots)
})

# Raw data split by cluster assignment ----

pred_plots <- lapply(PHENO, function (p) {
  res_plots <- lapply(STRATA, function (s) {
    if (!is.na(best_flex[[p]][[s]])) {
      cdat <- raw_dat[[p]][[s]]
      cdat$cluster <- as.factor(clusters(best_flex[[p]][[s]]))
      
      # Calculate mean and SE in each 5-year interval within each cluster
      age_bin_cuts <- seq(20, 80, by = 5)
      cdat$age_bin <- cut(cdat$age_event, 
                          age_bin_cuts, include.lowest = T)
      plot_df <- cdat %>% group_by(cluster, age_bin) %>% 
        summarise(count = n(),
                  mean_value = mean(value),
                  se_value = sd(value)/sqrt(count))
      
      # Plot 
      p <- ggplot(plot_df, aes(x = age_bin, y = mean_value,
                               group = cluster, 
                               color = cluster, fill = cluster)) +
        geom_point() +
        geom_path() +
        geom_ribbon(aes(ymin = mean_value - se_value, 
                        ymax = mean_value + se_value),
                    alpha = 0.2) +
        scale_fill_brewer(palette = "Set1") +
        scale_color_brewer(palette = "Set1") +
        labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
             title = s)
    } else { p <- NA }
    return (p)
  })
  pdf(paste0("plots/mixture_models/best_model_group_means_", p, ".pdf"), 
      onefile = T)
  print(res_plots)
  dev.off()
  return (res_plots)
})
