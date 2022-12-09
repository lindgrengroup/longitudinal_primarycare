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

# Report variance explained by each of the covariates used in baseline model ----

R2_table <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- raw_slopes[[p]][[s]]
    age_model <- lm(raw_slope ~ baseline_age + age_sq,
                    data = df)
    R2_age <- summary(age_model)$r.squared

    PC_model <- lm(formula(paste0("raw_slope ~ ", paste(PCs, collapse = " + "))),
                   data = df)
    R2_PCs <- summary(PC_model)$r.squared

    if (grepl("sexcomb", s)) {
      # Sex-combined analysis with sex as a covariate
      sex_model <- lm(raw_slope ~ sex, data = df)
      R2_sex <- summary(sex_model)$r.squared
    } else {
      R2_sex <- NA
    }
    res <- data.frame(adiposity = p,
                      ancestry = strsplit(s, "_")[[1]][1],
                      sex = strsplit(s, "_")[[1]][2],
                      R2_age = R2_age, R2_PCs = R2_PCs, R2_sex = R2_sex)
    return (res)
  })
  return (res)
})
R2_table <- bind_rows(R2_table)

write.table(R2_table, "results/adjusted_slopes/baseline_R2.txt",
            sep = "\t", row.names = F, quote = F)

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
# Don't adjust for height for BMI or WHR slopes 
formula_strings$m4bmi <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_BMI + FUyrs + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4bmisex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_BMI + FUyrs + ", 
         paste(PCs, collapse = " + "))

# Model 4_alt_1: adjust for baseline adipo-trait (not BMI)
formula_strings$m4a1 <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_trait + height + FUyrs + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a1sex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_trait + height + FUyrs + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a1bmi <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_trait + FUyrs + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a1bmisex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_trait + FUyrs + ", 
         paste(PCs, collapse = " + "))

# Model 4_alt_2: adjust for # FU measures as well as FU yrs
formula_strings$m4a2 <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_BMI + height + FUyrs + FU_n + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a2sex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_BMI + height + FUyrs + FU_n + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a2bmi <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_BMI + FUyrs + FU_n + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a2bmisex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_BMI + FUyrs + FU_n + ", 
         paste(PCs, collapse = " + "))

# Model 4_alt_3: adjust for baseline trait, FU_n as well as FU yrs
formula_strings$m4a3 <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_trait + height + FUyrs + FU_n + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a3sex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_trait + height + FUyrs + FU_n + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a3bmi <- 
  paste0("raw_slope ~ baseline_age + age_sq + baseline_trait + FUyrs + FU_n + ", 
         paste(PCs, collapse = " + "))
formula_strings$m4a3bmisex <- 
  paste0("raw_slope ~ baseline_age + age_sq + sex + baseline_trait + FUyrs + FU_n + ", 
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
    
    # Model 3 - don't run for BMI or WHR as it adjusts for height
    if (p == "BMI" | p == "WHR") {
      m3 <- NA 
    } else {
      if (grepl("sexcomb", s)) {
        m3 <- lm(formula(formula_strings$m3sex), df)
      } else m3 <- lm(formula(formula_strings$m3), df)
    }
    
    # Model 4
    if (p == "BMI" | p == "WHR") {
      if (grepl("sexcomb", s)) {
        m4 <- lm(formula(formula_strings$m4bmisex), df)
        m4a1 <- lm(formula(formula_strings$m4a1bmisex), df)
        m4a2 <- lm(formula(formula_strings$m4a2bmisex), df)
        m4a3 <- lm(formula(formula_strings$m4a3bmisex), df)
      } else {
        m4 <- lm(formula(formula_strings$m4bmi), df)
        m4a1 <- lm(formula(formula_strings$m4a1bmi), df)
        m4a2 <- lm(formula(formula_strings$m4a2bmi), df)
        m4a3 <- lm(formula(formula_strings$m4a3bmi), df)
      } 
    } else {
      if (grepl("sexcomb", s)) {
        m4 <- lm(formula(formula_strings$m4sex), df)
        m4a1 <- lm(formula(formula_strings$m4a1sex), df)
        m4a2 <- lm(formula(formula_strings$m4a2sex), df)
        m4a3 <- lm(formula(formula_strings$m4a3sex), df)
      } else {
        m4 <- lm(formula(formula_strings$m4), df)
        m4a1 <- lm(formula(formula_strings$m4a1), df)
        m4a2 <- lm(formula(formula_strings$m4a2), df)
        m4a3 <- lm(formula(formula_strings$m4a3), df)
      } 
    }
    
    return (list(m1 = m1, m2 = m2, m3 = m3, 
                 m4 = m4, m4a1 = m4a1, m4a2 = m4a2, m4a3 = m4a3))
  })
  names(res) <- STRATA
  return (res)
})
names(models) <- PHENOTYPES

saveRDS(models, "/well/lindgren/UKBIOBANK/samvida/adiposity/adj_slope_models.rds")

# Compare nested models with ANOVA ----

# Compare M1-M2-M3-M4
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

# Compare non-nested models with AIC ----

# Compare M4-alt1-alt2-alt3
lapply(PHENOTYPES, function (p) {
  lapply(STRATA, function (s) {
    ms <- models[[p]][[s]]
    sink(paste0("log_files/adj_slopes_model_AIC_", p, ".txt"), append = T)
    cat(paste0("Strata: ", s, " ", "\n"))
    print(AIC(ms$m4, ms$m4a1, ms$m4a2, ms$m4a3))
    cat("\n")
    sink()
  })
})


