# Author: Samvida S. Venkatesh
# Date: 14/12/20

library(tidyverse)
library(flexmix)
library(RColorBrewer)
theme_set(theme_bw())

SEED <- 1412
set.seed(SEED)

# Read cleaned data ----

# BMI <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/cleaned_all_BMI.txt",
#                   sep = "\t", header = T, stringsAsFactors = F)

BMI <- read.table("test_cleaned_BMI.txt",
                  sep = "\t", header = T, stringsAsFactors = F)

BMI$eid <- as.factor(as.character(BMI$eid))
BMI <- BMI[!is.na(BMI$BMI), ]
BMI <- BMI %>% group_by(eid) %>% arrange(age_years, .by_group = T)

# Rank-based inverse normal transformation of BMI values
BMImodel <- lm(BMI ~ 1, data = BMI, na.action = na.exclude)
BMI$BMI_RINT <- residuals(BMImodel)
BMI$BMI_RINT <- qnorm((rank(BMI$BMI_RINT, na.last = "keep") - 0.5) / 
                        sum(!is.na(BMI$BMI_RINT)))

# Group by number of measurements
BMI <- BMI %>% group_by(eid) %>% mutate(n_obs = n())
breakpoints <- c(-Inf, 5, 10, 15, Inf)
names <- c("2-5", "6-10", "11-15", "16+")
BMI$nobs_bin <- cut(BMI$n_obs, breaks = breakpoints, labels = names)

BMI <- BMI %>% group_by(nobs_bin) %>% group_split()

# Run baseline models (single-group) ----

g1 <- BMI[[1]]
g1 <- subset(g1, g1$sex == "F")

# Fit a series of models (no growth, linear, quadratic, )
baseline_nogrowth <- flexmix(.~ .|eid, k = 1, 
                              model = FLXMRlmm(BMI_RINT ~ 1,
                                               random = ~ 1),
                        data = g1)

baseline_linear <- flexmix(.~ .|eid, k = 1, 
                           model = FLXMRlmm(BMI_RINT ~ age_years,
                                            random = ~ 1),
                           data = g1)

baseline_quadratic <- flexmix(.~ .|eid, k = 1, 
                           model = FLXMRlmm(BMI_RINT ~ poly(age_years, 2),
                                            random = ~ 1),
                           data = g1)

baseline_cubic <- flexmix(.~ .|eid, k = 1, 
                         model = FLXMRlmm(BMI_RINT ~ poly(age_years, 3),
                                          random = ~ 1),
                         data = g1)

baseline_spline <- flexmix(.~ .|eid, k = 1, 
                           model = FLXMRlmm(BMI_RINT ~ age_years,
                                            random = ~ 1, 
                                            lm.fit = "smooth.spline"),
                           data = g1)


ids <- sample(unique(g1$eid), 100)
g1plot <- g1[g1$eid %in% ids, ]

fits <- data.frame(age_years = g1$age_years,
                   nogrowth = fitted(baseline_nogrowth),
                   linear = fitted(baseline_linear),
                   quadratic = fitted(baseline_quadratic),
                   cubic = fitted(baseline_cubic),
                   spline = fitted(baseline_spline))
colnames(fits) <- c("age_years", "nogrowth", "linear", "quadratic",
                    "cubic", "spline")

ggplot(g1plot) +
  geom_point(aes(x = age_years, y = BMI_RINT, group = eid), col = "#D3D3D3") +
  geom_line(aes(x = age_years, y = BMI_RINT, group = eid), col = "#D3D3D3") +
  geom_line(data = fits, 
            aes(x = age_years, y = nogrowth), col = "black") +
  geom_line(data = fits, 
            aes(x = age_years, y = linear), col = "#1B9E77") +
  geom_line(data = fits, 
            aes(x = age_years, y = quadratic), col = "#D95F02") +
  geom_line(data = fits, 
            aes(x = age_years, y = cubic), col = "#7570B3") +
  geom_line(data = fits, 
            aes(x = age_years, y = cubic), col = "#7570B3") +
  labs(x = "Age (years)", y = "BMI (R-INT)")

gmm1 <- stepFlexmix(.~ .|eid, k = 1:10, nrep = 100, 
                    model = FLXMRlmm(BMI_RINT ~ age_years, 
                                     random = ~ 1 + age_years,
                       varFix = c(Random=F, Residual=T)), 
                    data = g1)

gmm2_2 <- getModel(gmm1, which =2)

