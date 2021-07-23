# Author: Samvida S. Venkatesh
# Date: 19/04/2021

library(tidyverse)
theme_set(theme_bw())
library(flexmix)

set.seed(190421)

# Simulate data ----

# Simulate data from two normal distributions with different means, and different
# number of samples from each

# Set 1
m1 <- 0
m2 <- 50
sd1 <- 5
sd2 <- 5
N1 <- 100
N2 <- 10

# Set 2 - more challenging
m1 <- 20
m2 <- 40
sd1 <- 5
sd2 <- 10
N1 <- 100
N2 <- 100

# Set 3 - most challenging
m1 <- 21
m2 <- 22
sd1 <- 5
sd2 <- 10
N1 <- 100
N2 <- 100

a <- rnorm(n = N1, mean = m1, sd = sd1)
b <- rnorm(n = N2, mean = m2, sd = sd2)

# Mix the data by creating a single distribution
x <- c(a, b)
# Original class of data
class <- c(rep("a", N1), rep("b", N2))
# Bind to dataframe
dat <- data.frame(x = x, class = class)

# Plot simulated data
ggplot(dat, aes(x = x)) +
  geom_histogram(aes(x, ..density..), binwidth = 1, 
                 fill = "white", colour = "black") +
  geom_vline(xintercept = m1, colour = "red", size = 2) +
  geom_vline(xintercept = m2, colour = "blue", size = 2) 

# Model ----

mod1 <- FLXMRglm(family = "gaussian")
mod2 <- FLXMRglm(family = "gaussian")
flexfit <- flexmix(x ~ 1, data = dat, k = 2,
                   model = list(mod1, mod2))

# Look at cluster assignments
table(clusters(flexfit), dat$class)

# Look at parameters
c1 <- parameters(flexfit, component = 1)[[1]]
c2 <- parameters(flexfit, component = 2)[[1]]

# Plot mixture models on top of the real data
plot_mix_comps <- function (x, mu, sig, lam) {
  lam * dnorm(x, mu, sig)
}
# Scale height by number of points assigned to each component
lam <- table(clusters(flexfit))
ggplot(data = dat) +
  geom_histogram(aes(x, ..density..), binwidth = 1, 
                 fill = "white", colour = "black") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mu = c1[1], sig = c1[2], lam = lam[1]/sum(lam)),
                colour = "red", lwd = 1.5) + 
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mu = c2[1], sig = c2[2], lam = lam[2]/sum(lam)),
                colour = "blue", lwd = 1.5) + 
  ylab("density")

# Small subset of adiposity data for tutorial ----

# BMI, chinese F
dat <- readRDS("BMI_dat_subset.rds") %>% ungroup()

# Plot raw data
ggplot(dat, aes(x = age_event, y = value, group = eid)) +
  geom_point(colour = "#d3d3d3") +
  geom_line(colour = "#d3d3d3")

# Increasing class latent-class growth analysis models ----

# Model parameters: GLM, BMI ~ age with fixed parameters:
# variance in each component constrained to be equal, 
# i.e. single residual term that reflects the overall unexplained variance in data
# Can use FLXMRglm if no such constraint required
mod1 <- FLXMRglmfix(value ~ age_event, varFix = T)
# Run the model repeatedly for different number of components
# Parameters: use model defined above, nested within each individual
# 8 classes, nrep to keep solution with maximum likelihood
lcgaMix <- stepFlexmix(.~ .|eid, k = 1:4, nrep = 50, 
                       model = mod1, 
                       data = dat, control = list(iter.max = 500, minprior = 0))

# Inspect
lcgaMix

# Look at a specific class model
lcga4 <- getModel(lcgaMix, which = 4)
summary(lcga4)
parameters(lcga4)

# Plot
lcga4_par <- data.frame(parameters(lcga4))
age_range <- seq(min(dat$age_event), max(dat$age_event), length.out = 100)
lcga4_plot <- data.frame(age = age_range,
                         comp1 = lcga4_par$Comp.1[1] + 
                           age_range*lcga4_par$Comp.1[2],
                         comp2 = lcga4_par$Comp.2[1] + 
                           age_range*lcga4_par$Comp.2[2],
                         comp3 = lcga4_par$Comp.3[1] + 
                           age_range*lcga4_par$Comp.3[2],
                         comp4 = lcga4_par$Comp.4[1] + 
                           age_range*lcga4_par$Comp.4[2])
# Variance is fixed so we can just use one sigma
group_sig <- lcga4_par$Comp.1[3]
lcga4_plot <- pivot_longer(lcga4_plot, cols = starts_with("comp"),
                           names_to = "component", values_to = "value") %>%
  mutate(lci = value - 1.96*group_sig,
         uci = value + 1.96*group_sig)

ggplot(lcga4_plot, aes(x = age, y = value, 
                       colour = component, fill = component)) +
  geom_line() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  labs(x = "age", y = "BMI")

# Baseline model (single group) ----

# Linear mixed models with parameters for:
# fixed effect of age on adiposity trait, random slopes and intercepts 
# for each individual
mod1 <- FLXMRlmm(value ~ age_event, 
                 random = ~ 1 + age_event,
                 varFix = c(Random = F, Residual = F))
# Smoothing spline 
mod2 <- FLXMRlmm(value ~ 0 + age_event, 
                 random = ~ 1 + age_event, 
                 lm.fit = "smooth.spline",
                 varFix = c(Random = F, Residual = F))

# Flexmix for baseline models
flex_bl1 <- flexmix(.~.|eid, k = 1, model = mod1, 
                    data = dat, control = list(iter.max = 1000, minprior = 0))
summary(flex_bl1)

flex_bl2 <- flexmix(.~.|eid, k = 1, model = mod2, 
                    data = dat, control = list(iter.max = 1000, minprior = 0))
summary(flex_bl2)

# Compare models with AIC and BIC
AIC(flex_bl1, flex_bl2)
BIC(flex_bl1, flex_bl2)

# Visually compare models by plotting their results
# Sample 100 individuals for background trajectories to plot 
raw_dat <- sample(unique(dat$eid), 100, replace = F)
raw_dat <- subset(dat, dat$eid %in% raw_dat)

# Build data frame of model predictions
mod_preds <- function (x) {
  pred1 <- predict(flex_bl1, newdata = x,
                  se.fit = T)
  group_sd1 <- sqrt(parameters(flex_bl1)["sigma2.Residual", 1])
  res1 <- x %>% mutate(model = "linear",
                    mean = pred1$Comp.1[, 1],
                    lci = mean - group_sd1,
                    uci = mean + group_sd1)
  pred2 <- predict(flex_bl2, newdata = x,
                   se.fit = T)
  group_sd2 <- sqrt(parameters(flex_bl2)["sigma2.Residual", 1])
  res2 <- x %>% mutate(model = "spline",
                       mean = pred2$Comp.1[, 1],
                       lci = mean - group_sd2,
                       uci = mean + group_sd2)
  res <- bind_rows(res1, res2)
  return (res)
}

pred_dat <- data.frame(age_event = seq(min(dat$age_event), 
                                       max(dat$age_event), length.out = 100),
                       eid = 1)
pred_dat <- mod_preds(pred_dat)

ggplot(raw_dat, aes(x = age_event)) +
  geom_point(aes(y = value, group = eid), 
             colour = "#d3d3d3") +
  geom_line(aes(y = value, group = eid), 
            colour = "#d3d3d3", linetype = "dashed") +
  geom_line(data = pred_dat, 
            aes(x = age_event, y = mean, colour = model)) +
  geom_ribbon(data = pred_dat, 
              aes(x = age_event, ymin = lci, ymax = uci,
                  colour = model, fill = model),
              alpha = 0.2) +
  scale_colour_brewer(palette = "Set1") + 
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "Age", y = "Adiposity trait")
  
# Increasing number of clusters, step-wise ----

flex_step <- stepFlexmix(.~.|eid, k = 1:10, 
                         model = mod2, 
                         data = dat, 
                         control = list(iter.max = 1000, minprior = 0))

# BIC suggests that a 1-class model is best 
plot(1:10, BIC(flex_step))
best_flex <- getModel(flex_step, "BIC")


