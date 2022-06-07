# Author: Samvida S. Venkatesh
# Date: 15/06/21

library(lme4)
library(splines)
library(tidyverse)
theme_set(theme_bw())
library(merTools)

set.seed(150621)

# Read files ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/visually_QCd_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/model_covariates.rds")
PCs <- paste0("PC", 1:21)

# Adjust data for covariates ----

to_adj <- merge(adiposity, covars, by = "eid")

covar_form <- paste0("value ~ baseline_age + age_sq + FUyrs + FU_n + ",
                     paste0(PCs, collapse = " + "),
                     " + (1 | data_provider)")
covar_mod <- lmer(as.formula(covar_form), data = to_adj, REML = F)

# Perform all analyses on the residuals

adiposity$adj_value <- residuals(covar_mod)

plot_adipo <- subset(adiposity, adiposity$eid %in% plot_ids)
pdf("random_eids_adj_dat.pdf")
ggplot(plot_adipo, aes(x = age_event, y = adj_value)) +
  facet_wrap(~eid, ncol = 3) +
  geom_point() +
  labs(x = "Age (years)", y = "Adjusted BMI")
dev.off()

# Linear mixed effects models ----

# Model formulas
mod_linear_fixed <- "adj_value ~ age_event"
mod_quad_fixed <- "adj_value ~ age_event + I(age_event^2)"
mod_spline_fixed <- "adj_value ~ ns(age_event, 3)"

mod_random_intercept <- "(1 | eid)"
mod_random_linear <- "(1 + age_event | eid)"
mod_random_quad <- "(1 + age_event + I(age_event^2) | eid)"
mod_random_spline <- "(1 + ns(age_event, 3) | eid)"

## Run models sequentially ----

# Linear slope, random intercept
mod1 <- lmer(as.formula(paste(mod_linear_fixed, 
                              mod_random_intercept, sep = " + ")), 
             data = adiposity, REML = F)
# Linear slope, random slope and intercept
mod2 <- lmer(as.formula(paste(mod_linear_fixed, 
                              mod_random_linear, sep = " + ")), 
             data = adiposity, REML = F)

# Linear + quadratic slope, random intercept
mod3 <- lmer(as.formula(paste(mod_quad_fixed, 
                              mod_random_intercept, sep = " + ")), 
             data = adiposity, REML = F)
# Linear + quadratic slope, random linear slope and intercept
mod4 <- lmer(as.formula(paste(mod_quad_fixed, 
                              mod_random_linear, sep = " + ")), 
             data = adiposity, REML = F)
# Linear + quadratic slope, random linear + quad slope and intercept
mod5 <- lmer(as.formula(paste(mod_quad_fixed, 
                              mod_random_quad, sep = " + ")), 
             data = adiposity, REML = F)

# Spline slope, random intercept
mod6 <- lmer(as.formula(paste(mod_spline_fixed, 
                              mod_random_intercept, sep = " + ")), 
             data = adiposity, REML = F)

# Spline slope, random intercept + spline slope
mod7 <- lmer(as.formula(paste(mod_spline_fixed, 
                              mod_random_spline, sep = " + ")), 
             data = adiposity, REML = F)

## Compare models ----

modelNames <- c("~age, 1|ID", "~age, 1+age|ID",
                "~age+age^2, 1|ID", "~age+age^2, 1+age|ID", 
                "~age+age^2, 1+age+age^2|ID",
                "~ns(age,3), 1|ID", "~ns(age,3), 1+ns(age,3)|ID")
compare_df <- data.frame(modelname = rep(modelNames, times = 2),
                         IC_value = c(AIC(mod1, mod2, mod3, mod4, 
                                          mod5, mod6, mod7)$AIC,
                                      BIC(mod1, mod2, mod3, mod4, 
                                          mod5, mod6, mod7)$BIC),
                         IC_name = rep(c("AIC", "BIC"), 
                                       each = length(modelNames)))

compare_plot <- ggplot(compare_df, aes(x = modelname, y = IC_value,
                                       colour = IC_name)) +
  geom_point() +
  scale_colour_brewer(palette = "Dark2") + 
  labs(x = "model", y = "Info criterion") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

pdf("model_criteria.pdf", onefile = T)
print(compare_plot)
dev.off()

best_mod <- mod5

## Plot residuals ----

# Residuals vs fitted values
resid_dat <- data.frame(fitted = fitted(best_mod),
                        resid = residuals(best_mod))
prf <- ggplot(resid_dat, aes(x = fitted, y = resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "fitted value", y = "residual")

pdf("model_residuals.pdf", onefile = T)
# QQ-plot
qqnorm(residuals(best_mod))
abline(a = 0, b = 1)
print(prf)
dev.off()

# Get individual-level parameters ----

indiv_dat <- as.data.frame(ranef(best_mod)$eid)

# Cluster with k-means
kdat <- data.frame(eid = rownames(indiv_dat))
kss <- rep(0, 10)
for (k in 1:10) {
  k_alg <- kmeans(indiv_dat, k, nstart = 25)
  kdat$cluster <- k_alg$cluster
  colnames(kdat)[ncol(kdat)] <- paste0("cluster", k)
  kss[k] <- k_alg$tot.withinss
}
pdf("kmeans_elbow.pdf")
plot(1:10, kss)
dev.off()

# Decide on n clusters
plot_clustered_dat <- adiposity
plot_clustered_dat$cluster <- kdat$cluster5[match(plot_clustered_dat$eid,
                                                  kdat$eid)]
plot_clustered_dat$age_bin <- cut(plot_clustered_dat$age_event,
                                  seq(20, 80, by = 5))
plot_clustered_dat <- plot_clustered_dat %>% 
  group_by(cluster, age_bin) %>% 
  summarise(mean = mean(value),
            count = n(),
            lwr = mean(value) - (1.96*sd(value)/sqrt(count)),
            upr = mean(value) + (1.96*sd(value)/sqrt(count)))
plot_clustered_dat$cluster <- as.factor(plot_clustered_dat$cluster)

kplot <- ggplot(plot_clustered_dat, aes(x = age_bin, y = mean,
                                        colour = cluster, fill = cluster,
                                        group = cluster)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  scale_colour_brewer(palette = "Dark2", guide = F) +
  scale_fill_brewer(palette = "Dark2", guide = F) +
  labs(x = "Age bin (years)", y = "Mean (95% CI) BMI")

pdf("cluster_means.pdf", onefile = T)
print(kplot)
dev.off()

## Interpret regression coefficients ----

summary(best_mod)

# Plot population-level prediction
newdat <- data.frame(age_event = seq(20, 80, length.out = 100),
                     eid = 0)

pred_dat <- predictInterval(best_mod, newdata = newdat,
                            level = 0.95, n.sims = 1000,
                            stat = "median", type = "linear.prediction",
                            include.resid.var = TRUE)
pred_dat <- cbind(newdat, pred_dat)
p_FE <- ggplot(pred_dat, aes(x = age_event, y = fit)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2,
              colour = "grey") + 
  labs(x = "Age (years)", y = "Predicted adj. BMI")

# Plot individual-level predictions for a sample
newdat <- 
  data.frame(age_event = rep(seq(20, 80, length.out = 100), 
                             times = length(plot_ids)),
             eid = rep(plot_ids, each = 100))
pred_dat <- predictInterval(best_mod, newdata = newdat,
                            level = 0.95, n.sims = 1000,
                            stat = "median", type = "linear.prediction",
                            include.resid.var = TRUE)
pred_dat <- cbind(newdat, pred_dat)

p_RE <- ggplot(pred_dat, aes(x = age_event, y = fit)) +
  facet_wrap(~eid, ncol = 3) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = "grey") +
  labs(x = "Age (years)", y = "Adjusted BMI")

pdf("model_predictions.pdf", onefile = T)
print(p_FE)
print(p_RE)
dev.off()