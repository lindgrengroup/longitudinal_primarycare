# Author: Samvida S. Venkatesh
# Date: 15/06/21

library(lme4)
library(splines)
library(merTools)
library(tidyverse)
theme_set(theme_bw())

set.seed(150621)

PHENO <- "BMI"
SEX_STRATA <- "sex_comb"

# Read files ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/visually_QCd_adiposity.rds")
adiposity <- adiposity[[PHENO]]

covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/model_covariates.rds")
covars <- covars[[PHENO]]

dat <- merge(adiposity, covars, by = "eid")
dat <- subset(dat, dat$ancestry == "white")

if (SEX_STRATA != "sex_comb") {
  dat <- subset(dat, dat$sex == SEX_STRATA)
}

# Log file
log_file <- paste0("model_log_", PHENO, "_", SEX_STRATA, ".txt")

# Plot basis functions for non-linearity ----

ns3 <- cbind.data.frame(dat$age_event, ns(dat$age_event, 3)) 
colnames(ns3) <- c("age", "ns1", "ns2", "ns3")
ns3 <- ns3 %>% pivot_longer(names_to = "ns", values_to = "ns_value", 
                                cols = -age) 

ns3_plot <- ggplot(ns3, aes(x = age, y = ns_value, 
                                colour = ns)) +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Age", y = "Basis value", 
       title = "Natural cubic spline, 3 df")

# Mixed effects models ----

## Model formulas ----

PCs <- paste0("PC", 1:21)
if (SEX_STRATA == "sex_comb") {
  COVAR_TERMS <- c("sex", "baseline_age", "age_sq", "data_provider",
                   "FUyrs", "FU_n", "PCs")
} else {
  COVAR_TERMS <- c("baseline_age", "age_sq", "data_provider", 
                   "FUyrs", "FU_n", "PCs")
}

FE_TERMS <- c("no_growth", "linear_slope", "quad_slope", "ns_slope")
FE_FORMS <- c("value ~ 1", 
              "value ~ age_event",
              "value ~ age_event + I(age_event^2)",
              "value ~ ns(age_event, 3)")
FE_FORMS <- paste0(FE_FORMS, " + ",
                   paste0(COVAR_TERMS[-length(COVAR_TERMS)], 
                          collapse = " + "), " + ",
                   paste0(PCs, collapse = " + "), " + ",
                   "(1 | data_provider)")
names(FE_FORMS) <- FE_TERMS

RE_TERMS <- c("intercept", "linear_slope", "quad_slope", "ns_slope")
RE_FORMS <- c("(1 | eid)", 
              "(1 + age_event | eid)", 
              "(1 + age_event + I(age_event^2) | eid)", 
              "(1 + ns(age_event, 3) | eid)")
names(RE_FORMS) <- RE_TERMS

## Determine RE structure ----

# Start with elaborate FE to determine RE structure

re_models <- lapply(RE_FORMS, function (rf) {
  return (lmer(as.formula(paste(FE_FORMS[["ns_slope"]], 
                                rf, sep = " + ")), 
               data = dat, REML = F))
})

# Compare models
re_comp <- with(re_models, 
                do.call(anova, lapply(names(re_models), as.name)))

sink(log_file, append = T)
cat("Determining best random effects structure: ", "\n")
print(re_comp)
sink()

# Get the model with minimum AIC for RE structure 
# Favour more complex models (which is why we use AIC)
best_RE <- rownames(re_comp)[which.min(re_comp$AIC)]

## Determine FE structure ----

fe_models <- lapply(FE_FORMS, function (ff) {
  return (lmer(as.formula(paste(ff, 
                                RE_FORMS[[best_RE]], sep = " + ")), 
               data = dat, REML = F))
})

# Compare models
fe_comp <- with(fe_models, 
                do.call(anova, lapply(names(fe_models), as.name)))

sink(log_file, append = T)
cat("Determining best fixed effects structure: ", "\n")
print(fe_comp)
sink()

# Get the model with minimum BIC for FE structure 
best_FE <- rownames(fe_comp)[which.min(fe_comp$BIC)]

best_mod <- fe_models[[best_FE]]

# Save model summary
sink(log_file, append = T)
cat("Best model parameters: ", "\n")
summary(best_mod)
sink()

# Save model
saveRDS(best_mod, paste0("best_model_", PHENO, "_", SEX_STRATA, ".rds"))

# Plot residuals ----

# Residuals vs fitted values
resid_dat <- data.frame(fitted = fitted(best_mod),
                        resid = residuals(best_mod))

prf <- ggplot(resid_dat, aes(x = fitted, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "fitted value", y = "residual")

# Get model predictions ----

plot_ids <- sample(covars$eid, 12, replace = F)
N_AGE_PRED <- 100

# Individual-level prediction
newdat <- subset(covars, covars$eid %in% plot_ids)
newdat <- do.call("rbind", replicate(N_AGE_PRED, newdat, 
                                     simplify = FALSE))
newdat$age_event <- rep(seq(20, 80, length.out = 100),
                        each = length(plot_ids))
newdat$data_provider <- 0

pred_indiv <- predict(best_mod, newdata = newdat,
                      allow.new.levels = T)
pred_dat <- cbind(newdat, pred_indiv)
colnames(pred_dat)[ncol(pred_dat)] <- "fit"
pred_dat$fit_type <- "individual"

# Population-level prediction with bootstrapped 95% confidence intervals

newdat_pop <- newdat
newdat_pop$eid <- 0
newdat_pop$data_provider <- 0
pred_pop <- predictInterval(best_mod, newdata = newdat_pop,
                            level = 0.95, n.sims = 1000,
                            stat = "median", type = "linear.prediction",
                            include.resid.var = TRUE,
                            # ignore variance from covariate terms in FE
                            ignore.fixed.terms = 
                              c(COVAR_TERMS[-length(COVAR_TERMS)], PCs))
pred_pop <- cbind(newdat, pred_pop)
pred_pop$fit_type <- "population"

# Combine individual and population level predictions to plot
pred_dat <- bind_rows(pred_dat, pred_pop)
raw_dat <- subset(dat, dat$eid %in% plot_ids)

pred_plot <- ggplot(pred_dat, aes(x = age_event, y = fit)) +
  facet_wrap(~eid, ncol = 3) +
  # Show original data underneath
  geom_point(data = raw_dat, aes(x = age_event, y = value),
             colour = "black") +
  geom_ribbon(aes(ymin = lwr, ymax = upr, 
                  colour = fit_type, fill = fit_type), alpha = 0.2) + 
  geom_line(aes(colour = fit_type)) +
  labs(x = "Age (years)", y = paste(PHENO, "(fitted)"))

# Get random effect terms ---- 

indiv_dat <- as.data.frame(ranef(best_mod)$eid)
indiv_dat$eid <- rownames(indiv_dat)
# Append to covariates file
raw_slopes <- merge(covars, indiv_dat, by = "eid")

# Save table
write.table(raw_slopes, 
            paste0("full_slopes_", PHENO, "_", SEX_STRATA, ".txt"),
            sep = "\t", col.names = T, quote = F)

# All plots ----

pdf(paste0("model_plots_", PHENO, "_", SEX_STRATA, ".pdf"), onefile = T)
# Natural cubic spline basis
print(ns3_plot)
# QQ-plot
qqnorm(residuals(best_mod))
abline(a = 0, b = 1)
# Residuals vs fitted values
print(prf)
# Model predictions
print(pred_plot)
dev.off()
