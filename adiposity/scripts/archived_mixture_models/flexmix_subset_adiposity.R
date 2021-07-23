# Author: Samvida S. Venkatesh
# Date: 21/06/2021

library(lme4)
library(splines)
library(flexmix)
library(tidyverse)
theme_set(theme_bw())

set.seed(210621)

# Read data ----

dat <- readRDS("test_adj_BMI_subset.rds")

# Increasing class models for data ----

# Best model as determined earlier
model_form <- FLXMRlmm(adj_value ~ age_event + I(age_event^2), 
                       random = ~ 1 + age_event, 
                       varFix = c(Random = F, Residual = F))

flex_step <- stepFlexmix(.~.|eid, k = 1:10, 
                         model = model_form, 
                         data = dat, 
                         control = list(iter.max = 1000, minprior = 0))

# BIC suggests that a 1-class model is best 
plot(1:10, BIC(flex_step))
best_flex <- getModel(flex_step, "BIC")

summary(best_flex)

