# Author: Samvida S. Venkatesh
# Date: 11/12/20

library(tidyverse)
library(dplyr)
theme_set(theme_bw())
library(RColorBrewer)

# Read and supplement data ----

BMI <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/BMI_primarycare_indiv_clean.txt",
                  sep = "\t", header = T)

# Number of BMI measurements
BMI <- BMI %>% group_by(eid) %>% mutate(n_obs = n())

# Group by number of measurements
breakpoints <- c(-Inf, 1, 5, 10, 15, Inf)
names <- c("1", "2-5", "6-10", "11-15", "16+")
BMI$nobs_bin <- cut(BMI$n_obs, breaks = breakpoints, labels = names)

# Randomly sample trajectories in each nobs_bin ----

pids <- BMI %>% distinct(eid, sex, nobs_bin) %>%
  group_by(sex, nobs_bin) %>%
  sample_n(size = min(n(), 20))

pdata <- BMI[BMI$eid %in% pids$eid, ]

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/trajectories/sample_by_nobs.pdf",
    onefile = T)

# Females
ggplot(subset(pdata, pdata$sex == "F"), 
       aes(x = age_years, y = primarycare_BMI, group = eid)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point(col = "#F8766D") +
  geom_line(col = "#F8766D", size = 0.7) +
  labs(x = "Age (years)", y = "BMI")

# Males
ggplot(subset(pdata, pdata$sex == "M"), 
       aes(x = age_years, y = primarycare_BMI, group = eid)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point(col = "#00BFC4") +
  geom_line(col = "#00BFC4", size = 0.7) +
  labs(x = "Age (years)", y = "BMI")

dev.off()

# Randomly sample trajectories in each nobs_bin and overall change in BMI ----

# Classify individuals as "overall gain", "overall loss", or "stable"
# based on whether the change from initial to final BMI is > or < 10%

change <- subset(BMI, BMI$n_obs > 1) %>% group_by(eid) %>%
  mutate(change_percent = ((last(primarycare_BMI) - first(primarycare_BMI)) / 
           first(primarycare_BMI)) * 100,
         change_class = ifelse(change_percent > 10, "gain (>10%)",
                               ifelse(change_percent < -10, "loss (<-10%)",
                                      "stable (+/-10%)")))

pids <- change %>% group_by(sex, nobs_bin, change_class) %>%
  sample_n(size = min(n(), 10))

pdata <- change[change$eid %in% pids$eid, ]

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/trajectories/sample_by_change.pdf",
    onefile = T)

# Females
ggplot(subset(pdata, pdata$sex == "F"), 
       aes(x = age_years, y = primarycare_BMI, group = eid, 
           col = change_class)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point() +
  geom_line(size = 0.7) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Age (years)", y = "BMI", col = "BMI change",
       title = "Females")

# Males
ggplot(subset(pdata, pdata$sex == "M"), 
       aes(x = age_years, y = primarycare_BMI, group = eid, 
           col = change_class)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point() +
  geom_line(size = 0.7) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Age (years)", y = "BMI", col = "BMI change",
       title = "Males")

dev.off()





