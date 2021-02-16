# Author: Samvida S. Venkatesh
# Date: 11/12/20

library(tidyverse)
library(dplyr)
library(eeptools)
theme_set(theme_bw())
library(RColorBrewer)

# Read primary care data ----

BMI <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/BMI_primarycare_indiv_clean.txt",
                  sep = "\t", header = T, stringsAsFactors = F)

# Add UKBIOBANK BMI data ----

pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

valid_ids <- unique(BMI$eid)
pheno <- subset(pheno, pheno$f.eid %in% valid_ids)

# Collect fields -
# Date of attending assessment centre:
# f53-0.0, f53-1.0, f53-2.0
# BMI:
# f21001-0.0, f21001-2.0, f21001-2.0
# Sex: f31-0.0
# Year and month of birth: f34-0.0, f52-0.0

UKBB_BMI <- pheno[, c("f.eid", "f.53.0.0", "f.53.1.0", "f.53.2.0",
                      "f.21001.0.0", "f.21001.1.0", "f.21001.2.0",
                      "f.31.0.0", "f.34.0.0", "f.52.0.0")]
colnames(UKBB_BMI) <- c("eid", "e1", "e2", "e3",
                        "v1", "v2", "v3",
                        "sex", "year_of_birth", "month_of_birth")
UKBB_BMI$sex <- ifelse(UKBB_BMI$sex == 0, "F", "M")
UKBB_BMI$dob <- as.Date(paste(UKBB_BMI$year_of_birth, 
                              UKBB_BMI$month_of_birth,
                              1, sep = "-"))

UKBB_BMI <- pivot_longer(UKBB_BMI, cols = 2:7,
                         names_to = c(".value", "tmp"),
                         names_pattern = "(.)(.)",
                         values_drop_na = T)

UKBB_BMI <- UKBB_BMI[, c("eid", "sex", "dob", "e", "v")]
colnames(UKBB_BMI) <- c("eid", "sex", "dob", "event_dt", "BMI")
UKBB_BMI <- UKBB_BMI[!is.na(UKBB_BMI$BMI), ]

# Calculate age (although we have the age recorded at assessment centre,
# the precision is to year, while we can calculate month)
UKBB_BMI$dob <- as.Date(UKBB_BMI$dob)
UKBB_BMI$event_dt <- as.Date(UKBB_BMI$event_dt)

UKBB_BMI$age_years <- age_calc(UKBB_BMI$dob, UKBB_BMI$event_dt, units = "years")

UKBB_BMI <- UKBB_BMI %>% group_by(eid) %>% mutate(mean_UKBB_BMI = mean(BMI))

# Merge with GP data

colnames(BMI) <- c("eid", "event_dt", "BMI", "data_provider",
                   "sex", "age_years", "mean_UKBB_BMI")
BMI$event_dt <- as.Date(BMI$event_dt)

BMI <- merge(BMI, 
             UKBB_BMI[, c("eid", "event_dt", "BMI", "sex",
                          "age_years", "mean_UKBB_BMI")],
             all = T)

write.table(BMI, "/well/lindgren/UKBIOBANK/samvida/BMI/cleaned_all_BMI.txt",
            sep = "\t", quote = F, row.names = F)

# Read cleaned and supplemented data ----

BMI <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/cleaned_all_BMI.txt",
                  sep = "\t", header = T, stringsAsFactors = F)

# Number of BMI measurements
BMI <- BMI %>% group_by(eid) %>% mutate(n_obs = n())

# Group by number of measurements
breakpoints <- c(-Inf, 5, 10, 15, Inf)
names <- c("2-5", "6-10", "11-15", "16+")
BMI$nobs_bin <- cut(BMI$n_obs, breaks = breakpoints, labels = names)

# Rank-based inverse normal transform the data
BMImodel <- lm(BMI ~ 1, data = BMI, na.action = na.exclude)
BMI$BMI_RINT <- residuals(BMImodel)
BMI$BMI_RINT <- qnorm((rank(BMI$BMI_RINT, na.last = "keep") - 0.5) / 
                        sum(!is.na(BMI$BMI_RINT)))

# Randomly sample trajectories in each nobs_bin ----

pids <- BMI %>% distinct(eid, sex, nobs_bin) %>%
  group_by(sex, nobs_bin) %>%
  sample_n(size = min(n(), 20))

pdata <- BMI[BMI$eid %in% pids$eid, ]

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/trajectories/sample_by_nobs.pdf",
    onefile = T)

# BMI

# Females
ggplot(subset(pdata, pdata$sex == "F"), 
       aes(x = age_years, y = BMI, group = eid)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point(col = "#F8766D") +
  geom_line(col = "#F8766D", size = 0.7) +
  labs(x = "Age (years)", y = "BMI")

# Males
ggplot(subset(pdata, pdata$sex == "M"), 
       aes(x = age_years, y = BMI, group = eid)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point(col = "#00BFC4") +
  geom_line(col = "#00BFC4", size = 0.7) +
  labs(x = "Age (years)", y = "BMI")

# BMI-INT

ggplot(subset(pdata, pdata$sex == "F"), 
       aes(x = age_years, y = BMI_RINT, group = eid)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point(col = "#F8766D") +
  geom_line(col = "#F8766D", size = 0.7) +
  labs(x = "Age (years)", y = "BMI-RINT")

# Males
ggplot(subset(pdata, pdata$sex == "M"), 
       aes(x = age_years, y = BMI_RINT, group = eid)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point(col = "#00BFC4") +
  geom_line(col = "#00BFC4", size = 0.7) +
  labs(x = "Age (years)", y = "BMI-RINT")


dev.off()

# Randomly sample trajectories in each nobs_bin and overall change in BMI ----

# Classify individuals as "overall gain", "overall loss", or "stable"
# based on whether the change from initial to final BMI is > or < 10%

change <- BMI %>% group_by(eid) %>%
  arrange(age_years, .by_group = T) %>%
  mutate(change_percent = ((last(BMI) - first(BMI)) / 
                             first(BMI)) * 100,
         change_class = ifelse(change_percent > 10, "gain (>10%)",
                               ifelse(change_percent < -10, "loss (<-10%)",
                                      "stable (+/-10%)")))

pids <- change %>% group_by(sex, nobs_bin, change_class) %>%
  sample_n(size = min(n(), 10))

pdata <- change[change$eid %in% pids$eid, ]

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/trajectories/sample_by_change.pdf",
    onefile = T)

# BMI

# Females
ggplot(subset(pdata, pdata$sex == "F"), 
       aes(x = age_years, y = BMI, group = eid, 
           col = change_class)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point() +
  geom_line(size = 0.7) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Age (years)", y = "BMI", col = "BMI change",
       title = "Females")

# Males
ggplot(subset(pdata, pdata$sex == "M"), 
       aes(x = age_years, y = BMI, group = eid, 
           col = change_class)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point() +
  geom_line(size = 0.7) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Age (years)", y = "BMI", col = "BMI change",
       title = "Males")

# BMI-RINT

# Females
ggplot(subset(pdata, pdata$sex == "F"), 
       aes(x = age_years, y = BMI_RINT, group = eid, 
           col = change_class)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point() +
  geom_line(size = 0.7) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Age (years)", y = "BMI-RINT", col = "BMI change",
       title = "Females")

# Males
ggplot(subset(pdata, pdata$sex == "M"), 
       aes(x = age_years, y = BMI_RINT, group = eid, 
           col = change_class)) +
  facet_wrap(~nobs_bin, nrow = 3) +
  geom_point() +
  geom_line(size = 0.7) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Age (years)", y = "BMI-RINT", col = "BMI change",
       title = "Males")

dev.off()
