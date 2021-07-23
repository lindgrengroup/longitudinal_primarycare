# Author: Samvida S. Venkatesh
# Date: 28/06/21

library(tidyverse)
theme_set(theme_bw())

PHENO <- "BMI"
SEX_STRATA <- "F"

# Read files ----

dat <- read.table(paste0("raw_slopes_", PHENO, "_", SEX_STRATA, ".txt"),
                  sep = "\t", header = T)

log_file <- paste0("adj_slopes_log_", PHENO, "_", SEX_STRATA, ".txt")

PCs <- paste0("PC", 1:21)
if (SEX_STRATA == "sex_comb") {
  COVARS <- c("sex", "baseline_age", "age_sq", "FUyrs", "FU_n", "PCs")
} else {
  COVARS <- c("baseline_age", "age_sq", "FUyrs", "FU_n", "PCs")
}

# Currently there are 32 covariate / other info columns
# This may change so make sure to check
SLOPE_TERMS <- colnames(dat)[32:ncol(dat)]

# Proportion of variance in slope terms explained by each covariate ----

## Calculate linear models ----

term_var_explained <- function (slope_term) {
  res_list <- lapply(COVARS, function (cov) {
    if (cov == "PCs") {
      mod_form <- paste0(slope_term, " ~ ", paste0(PCs, collapse = " + "))
    } else {
      mod_form <- paste0(slope_term, " ~ ", cov)
    }
    mod_R2 <- summary(lm(as.formula(mod_form), data = dat))$r.squared
    return (mod_R2)
  })
  names(res_list) <- COVARS
  res <- bind_cols(res_list)
}

var_explained <- lapply(SLOPE_TERMS, function (s) term_var_explained(s) )
var_explained <- bind_rows(var_explained)
var_explained$term <- SLOPE_TERMS

sink(log_file, append = T)
cat("Proportion of variance explained by:", "\n")
print(var_explained)
cat("\n")
sink()

## Plot proportion of variance explained ----

plot_df <- pivot_longer(var_explained, cols = -term,
                        names_to = "covariate",
                        values_to = "R2")

covar_plot <- ggplot(plot_df, aes(x = term, y = R2,
                                  fill = covariate)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set3") +
  coord_flip() 

pdf(paste0("covar_R2_", PHENO, "_", SEX_STRATA, ".pdf"))
print(covar_plot)
dev.off()

# Full covariate models ----

full_cov_models <- lapply(SLOPE_TERMS, function (s) {
  mod_form <- paste0(s, " ~ ", 
                     paste0(COVARS[which(COVARS != "PCs")], 
                            collapse = " + "), " + ",
                     paste0(PCs, collapse = " + "))
  mod <- lm(as.formula(mod_form), data = dat)
  
  # Log variance explained
  sink(log_file, append = T)
  cat(paste0("Proportion of variance in ", s, 
             " explained by all covariates:"),
      "\n")
  print(summary(mod)$r.squared)
  cat("\n")
  sink()
  
  # Save model
  return (mod)
}) 
names(full_cov_models) <- SLOPE_TERMS

## Adjust coefficient terms for covariates ----

adj_slopes <- lapply(SLOPE_TERMS, function (s) {
  return (full_cov_models[[s]]$residuals)
})
names(adj_slopes) <- paste0("adj_", SLOPE_TERMS)
adj_slopes <- bind_cols(adj_slopes)
adj_slopes$eid <- dat$eid

# Save table
write.table(adj_slopes, 
            paste0("adj_slopes_", PHENO, "_", SEX_STRATA, ".txt"),
            sep = "\t", col.names = T, quote = F)

