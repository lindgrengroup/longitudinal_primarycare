# Author: Samvida S. Venkatesh
# Date: 02/03/21
# ADAPTED FROM ADIPOSITY CHANGE CONSORTIUM SOP

library(tidyverse)

# Read stratified raw slope and covariate data ----

raw_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_raw_slopes.rds")

PHENOTYPES <- names(raw_slopes)
NPCs <- 21
PCs <- paste0("PC", 1:NPCs)

# Create sex-combined dataframe within each ancestry group ----

# Ensure EIDs are characters to enable binding rows
raw_slopes <- lapply(PHENOTYPES, function (p) {
  res <- lapply(raw_slopes[[p]], function (x) {
    x$eid <- as.character(x$eid)
    return (x)
  })
  return (res)
})
names(raw_slopes) <- PHENOTYPES

# Bind rows within the same ancestry
raw_slopes <- lapply(PHENOTYPES, function (p) {
  temp <- bind_rows(raw_slopes[[p]])
  temp <- split(temp, f = temp$ancestry)
  names(temp) <- lapply(temp, function (x) paste(unique(x$ancestry), "sexcomb", 
                                                 sep = "_") )
  temp <- c(raw_slopes[[p]], temp)
}) 
names(raw_slopes) <- PHENOTYPES 
STRATA <- names(raw_slopes[[1]])

# Save new slopes list with sex-combined dataframes
saveRDS(raw_slopes, "/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_raw_slopes.rds")

# # Report variance explained by each of the covariates used in baseline model ----
# 
# R2_table <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(STRATA, function (s) {
#     df <- raw_slopes[[p]][[s]]
#     age_model <- lm(raw_slope ~ baseline_age + age_sq, 
#                     data = df)
#     R2_age <- summary(age_model)$r.squared
#     
#     PC_model <- lm(formula(paste0("raw_slope ~ ", paste(PCs, collapse = " + "))), 
#                    data = df)
#     R2_PCs <- summary(PC_model)$r.squared
#     
#     if (grepl("sexcomb", s)) {
#       # Sex-combined analysis with sex as a covariate
#       sex_model <- lm(raw_slope ~ sex, data = df)
#       R2_sex <- summary(sex_model)$r.squared
#     } else {
#       R2_sex <- NA
#     }
#     res <- data.frame(adiposity = p, 
#                       ancestry = strsplit(s, "_")[[1]][1],
#                       sex = strsplit(s, "_")[[1]][2],
#                       R2_age = R2_age, R2_PCs = R2_PCs, R2_sex = R2_sex)
#     return (res)
#   })
#   return (res)
# })
# R2_table <- bind_rows(R2_table)
# 
# write.table(R2_table, "results/adjusted_slopes/baseline_R2.txt",
#             sep = "\t", row.names = F, quote = F)

# Build nested models of increasing complexity ----

formula_strings <- list()
# Model 1 covariates: baseline age and age-sq, genetic PCs, (sex)
formula_strings$m1 <- paste0("raw_slope ~ baseline_age + age_sq + ", 
                             paste(PCs, collapse = " + "))
formula_strings$m1sex <- paste0("raw_slope ~ baseline_age + age_sq + sex + ", 
                                paste(PCs, collapse = " + "))
# Model 2: add baseline BMI
formula_strings$m2 <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_BMI + ", 
         paste(PCs, collapse = " + "))
formula_strings$m2sex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_BMI + ", 
         paste(PCs, collapse = " + "))

# Model 3: add median height
formula_strings$m3 <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_BMI + height + ", 
         paste(PCs, collapse = " + "))
formula_strings$m3sex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_BMI + height + ", 
         paste(PCs, collapse = " + "))

# Model 4: add FU years
formula_strings$m4 <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_BMI + height + FUyrs + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4sex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_BMI + height + FUyrs + ", 
         paste(PCs, collapse = " + "))
# Don't adjust for height for BMI slopes 
formula_strings$m4bmi <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_BMI + FUyrs + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4bmisex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_BMI + FUyrs + ", 
         paste(PCs, collapse = " + "))

# Run models ----

models <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- raw_slopes[[p]][[s]]
    
    # Model 1
    if (grepl("sexcomb", s)) {
      m1 <- lm(formula(formula_strings$m1sex), df)
    } else m1 <- lm(formula(formula_strings$m1), df)
    
    # Model 2
    if (grepl("sexcomb", s)) {
      m2 <- lm(formula(formula_strings$m2sex), df)
    } else m2 <- lm(formula(formula_strings$m2), df)
    
    # Model 3 - don't run for BMI as it adjusts for height
    if (p == "BMI") {
      m3 <- NA 
    } else {
      if (grepl("sexcomb", s)) {
        m3 <- lm(formula(formula_strings$m3sex), df)
      } else m3 <- lm(formula(formula_strings$m3), df)
    }
    
    # Model 4
    if (p == "BMI") {
      if (grepl("sexcomb", s)) {
        m4 <- lm(formula(formula_strings$m4bmisex), df)
      } else m4 <- lm(formula(formula_strings$m4bmi), df)
    } else {
      if (grepl("sexcomb", s)) {
        m4 <- lm(formula(formula_strings$m4sex), df)
      } else m4 <- lm(formula(formula_strings$m4), df)
    }
    
    return (list(m1 = m1, m2 = m2, m3 = m3, m4 = m4))
  })
  names(res) <- STRATA
  return (res)
})
names(models) <- PHENOTYPES

saveRDS(models, "/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_adj_slope_models.rds")

# Compare models with ANOVA ----

lapply(PHENOTYPES, function (p) {
  lapply(STRATA, function (s) {
    ms <- models[[p]][[s]]
    sink(paste0("log_files/adj_slopes_model_anova_", p, ".txt"), append = T)
    cat(paste0("Strata: ", s, " ", "\n"))
    if (p == "BMI") {
      print(anova(ms$m1, ms$m2, ms$m4))
    } else {
      print(anova(ms$m1, ms$m2, ms$m3, ms$m4))
    }
    cat("\n")
    sink()
  })
})

