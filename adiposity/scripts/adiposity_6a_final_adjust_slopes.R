# Author: Samvida S. Venkatesh
# Date: 02/03/21
# ADAPTED FROM ADIPOSITY CHANGE CONSORTIUM SOP

library(tidyverse)

# Read stratified raw slope and covariate data ----

raw_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/raw_slopes_and_covars.rds")

PHENOTYPES <- names(raw_slopes)
STRATA <- names(raw_slopes[[1]])
NPCs <- 21
PCs <- paste0("PC", 1:NPCs)

# Build final models ----

formula_strings <- list()
# Model covariates for BMI and WHR: 
# baseline age and age-sq, baseline trait, 
# # FU years, # FU measures,
# genetic PCs, (sex)
formula_strings$mA <- paste0("raw_slope ~ baseline_age + age_sq + baseline_trait + FUyrs + FU_n + ", 
                             paste(PCs, collapse = " + "))
formula_strings$mAsex <- paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_trait + FUyrs + FU_n + ", 
                                paste(PCs, collapse = " + "))

# Model covariates for weight and waist circumference: 
# baseline age and age-sq, baseline trait, median adult height,
# # FU years, # FU measures, 
# genetic PCs, (sex)
formula_strings$mB <- paste0("raw_slope ~ baseline_age + age_sq + baseline_trait + height + FUyrs + FU_n + ", 
                             paste(PCs, collapse = " + "))
formula_strings$mBsex <- paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_trait + height + FUyrs + FU_n + ", 
                                paste(PCs, collapse = " + "))

# Run models ----

models <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- raw_slopes[[p]][[s]]
    
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

saveRDS(models, "/well/lindgren/UKBIOBANK/samvida/adiposity/final_adj_model.rds")
