# Author: Samvida S. Venkatesh
# Date: 13/01/2021

library(tidyverse)
theme_set(theme_bw())

# Read data ----

pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

# Remove individuals who have withdrawn consent
withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20200820.csv", 
                        header = F)

pheno <- subset(pheno, !pheno$f.eid %in% withdrawn$V1)

case_control_ids <- readRDS("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH_FSH_case_control_ids.rds")
# names(case_control_ids) <- c("case_eids", "control_eids_s1", 
#                              "control_eids_s2", "control_eids_s3")

# Keep columns of interest in UKB data ----

# 1. Sex: f.31.0.0 (0 - F, 1 - M)
# 2. Birth year: f.34.0.0
# 3. BMI: f.21001.0.0, f.21001.1.0, f.21001.2.0
# 4a. Past smoker: f.1249.0.0, f.1249.1.0, f.1249.2.0 (0 - No, 1 & 2 - Yes)
# 4b. Current smoker: f.1239.0.0, f.1239.1.0, f.1239.2.0 (1 & 2 - Yes, 3 & 4 - No)
# 5. Alcohol intake: f.1558.0.0, f.1558.1.0, f.1558.2.0
# (1, 2, & 3 - more than once a week, 4 & 5 - less than once a week, 6 - never)
# 6. Ever had stillbirth/miscarriage/spontaneous termination: 
# f.2774.0.0, f.2774.1.0, f.2774.2.0 (0 - No, 1 - Yes)
# 7. Number of:
#   a. stillbirths - f.3829.0.0, f.3829.1.0, f.3829.2.0
#   b. spontaneous miscarriages - f.3839.0.0, f.3839.1.0, f.3839.2.0
#   c. pregnancy terminations - f.3849.0.0, f.3849.1.0, f.3849.2.0
# 8. Number of:
#   a. children fathered - f.2405.0.0, f.2405.1.0, f.2405.2.0
#   b. live births - f.2734.0.0, f.2734.1.0, f.2734.2.0
# 9. Ever taken oral contraceptive: f.2784.0.0, f.2784.1.0, f.2784.2.0 
# (0 - No, 1 - Yes)
# 10. Ever used hormone replacement therapy - f.2814.0.0, f.2814.1.0, f.2814.2.0 
# (0 - No, 1 - Yes)
# 11. Recalled age at:
#   a. menarche - f.2714.0.0, f.2714.1.0, f.2714.2.0 
#   b. menopause - f.3581.0.0, f.3581.1.0, f.3581.2.0
# 12. Ever had hysterectomy: f.3581.0.0, f.3581.1.0, f.3581.2.0

keep_ids <- unlist(case_control_ids)
cleaned <- subset(pheno, pheno$f.eid %in% keep_ids)

cleaned <- cleaned[, c("f.eid",
                       "f.31.0.0",
                       "f.34.0.0", 
                       "f.21001.0.0", "f.21001.1.0", "f.21001.2.0",
                       "f.1249.0.0", "f.1249.1.0", "f.1249.2.0",
                       "f.1239.0.0", "f.1239.1.0", "f.1239.2.0",
                       "f.1558.0.0", "f.1558.1.0", "f.1558.2.0",
                       "f.2774.0.0", "f.2774.1.0", "f.2774.2.0",
                       "f.3829.0.0", "f.3829.1.0", "f.3829.2.0",
                       "f.3839.0.0", "f.3839.1.0", "f.3839.2.0",
                       "f.3849.0.0", "f.3849.1.0", "f.3849.2.0",
                       "f.2405.0.0", "f.2405.1.0", "f.2405.2.0",
                       "f.2734.0.0", "f.2734.1.0", "f.2734.2.0",
                       "f.2784.0.0", "f.2784.1.0", "f.2784.2.0",
                       "f.2814.0.0", "f.2814.1.0", "f.2814.2.0",
                       "f.2714.0.0", "f.2714.1.0", "f.2714.2.0",
                       "f.3581.0.0", "f.3581.1.0", "f.3581.2.0",
                       "f.3591.0.0", "f.3591.1.0", "f.3591.2.0")]
colnames(cleaned) <- c("eid", 
                       "sex",
                       "year_of_birth",
                       "BMI_1", "BMI_2", "BMI_3",
                       "past_smoker_1", "past_smoker_2", "past_smoker_3",
                       "current_smoker_1", "current_smoker_2", "current_smoker_3",
                       "alc_intake_1", "alc_intake_2", "alc_intake_3",
                       "SSMT_1", "SSMT_2", "SSMT_3",
                       "stillbirths_1", "stillbirths_2", "stillbirths_3", 
                       "miscarr_1", "miscarr_2", "miscarr_3",
                       "term_1", "term_2", "term_3", 
                       "children_fathered_1", "children_fathered_2", "children_fathered_3", 
                       "live_births_1", "live_births_2", "live_births_3", 
                       "oral_contraceptive_1", "oral_contraceptive_2", "oral_contraceptive_3", 
                       "HRT_1", "HRT_2", "HRT_3",
                       "age_menarche_1", "age_menarche_2", "age_menarche_3", 
                       "age_menopause_1", "age_menopause_2", "age_menopause_3", 
                       "hysterectomy_1", "hysterectomy_2", "hysterectomy_3")

# Create columns of interest by combining the columns above

dat <- data.frame(eid = cleaned$eid)

dat$sex <- ifelse(cleaned$sex == 0, "F", "M")

dat$age_2020 <- 2020 - cleaned$year_of_birth

dat$BMI <- rowMeans(cleaned %>% select(starts_with("BMI")), na.rm = T)

dat$past_smoker <- apply(cleaned %>% select(starts_with("past_smoker")),
                         1, function (r) {
                           ifelse(any(r %in% c(1, 2), na.rm = T), T,
                                  ifelse(any(r %in% c(3, 4), na.rm = T), F, 
                                         NA))
                         })

dat$current_smoker <- apply(cleaned %>% select(starts_with("current_smoker")),
                            1, function (r) {
                              ifelse(any(r %in% c(1, 2), na.rm = T), T,
                                     ifelse(any(r == 0, na.rm = T), F, 
                                            NA))
                            })

dat$alc_intake <- apply(cleaned %>% select(starts_with("alc_intake")),
                        1, function (r) {
                          ifelse(any(r %in% c(1, 2, 3), na.rm = T), 
                                 "frequent",
                                 ifelse(any(r %in% c(4, 5), na.rm = T), 
                                        "occasional", 
                                        ifelse(any(r == 6, na.rm = T), 
                                               "never", NA)))
                        })

dat$SSMT <- apply(cleaned %>% select(starts_with("SSMT")),
                  1, function (r) {
                    ifelse(any(r == 1, na.rm = T), T,
                           ifelse(any(r == 0, na.rm = T), F, 
                                  NA))
                  })

dat$n_stillbirths <- apply(cleaned %>% select(starts_with("stillbirths")),
                           1, function (r) { 
                             ifelse(all(is.na(r)), NA, 
                                    ifelse(any(r < 0, na.rm = T), NA,
                                           median(r, na.rm = T)))
                           })

dat$n_miscarr <- apply(cleaned %>% select(starts_with("miscarr")),
                       1, function (r) { 
                         ifelse(all(is.na(r)), NA, 
                                ifelse(any(r < 0, na.rm = T), NA,
                                       median(r, na.rm = T)))
                       })

dat$n_terminations <- apply(cleaned %>% select(starts_with("term")),
                            1, function (r) { 
                              ifelse(all(is.na(r)), NA, 
                                     ifelse(any(r < 0, na.rm = T), NA,
                                            median(r, na.rm = T)))
                            })

dat$n_children_fathered <- apply(cleaned %>% select(starts_with("children")),
                                 1, function (r) { 
                                   ifelse(all(is.na(r)), NA, 
                                          ifelse(any(r < 0, na.rm = T), NA,
                                                 median(r, na.rm = T)))
                                 })

dat$n_live_births <- apply(cleaned %>% select(starts_with("live")),
                           1, function (r) { 
                             ifelse(all(is.na(r)), NA, 
                                    ifelse(any(r < 0, na.rm = T), NA,
                                           median(r, na.rm = T)))
                           })

dat$ever_oral_contraceptive <- apply(cleaned %>% 
                                       select(starts_with("oral_contraceptive")),
                                     1, function (r) {
                                       ifelse(any(r == 1, na.rm = T), T,
                                              ifelse(any(r == 0, na.rm = T), F, 
                                                     NA))
                                     })

dat$ever_HRT <- apply(cleaned %>% select(starts_with("HRT")),
                      1, function (r) {
                        ifelse(any(r == 1, na.rm = T), T,
                               ifelse(any(r == 0, na.rm = T), F, 
                                      NA))
                      })

dat$age_menarche <- apply(cleaned %>% select(starts_with("age_menarche")),
                          1, function (r) { 
                            ifelse(all(is.na(r)), NA, 
                                   ifelse(any(r < 0, na.rm = T), NA,
                                          median(r, na.rm = T)))
                          })

dat$age_menopause <- apply(cleaned %>% select(starts_with("age_menopause")),
                           1, function (r) { 
                             ifelse(all(is.na(r)), NA, 
                                    ifelse(any(r < 0, na.rm = T), NA,
                                           median(r, na.rm = T)))
                           })

dat$hysterectomy <- apply(cleaned %>% select(starts_with("hysterectomy")),
                          1, function (r) {
                            ifelse(any(r == 1, na.rm = T), T,
                                   ifelse(any(r == 0, na.rm = T), F, 
                                          NA))
                          })

# Calculate demographics on data ----

dat$case <- dat$eid %in% case_control_ids$case_eids
dat$control_s1 <- dat$eid %in% case_control_ids$control_eids_s1
dat$control_s2 <- dat$eid %in% case_control_ids$control_eids_s2
dat$control_s3 <- dat$eid %in% case_control_ids$control_eids_s3

write.table(dat, "summarised_LH_FSH_case_control_UKB_pheno.txt", sep = "\t",
            quote = F, row.names = F)

# Repeat this for each group (cases, all three control sets)

d <- subset(dat, dat$control_s3) 

# Sex
d %>% group_by(sex) %>% count()

# Group by sex before calculating other characteristics
# Quantitative
quant_summary <- d %>% group_by(sex) %>% 
  summarise(age_2020_mean = mean(age_2020, na.rm = T),
            age_2020_sd = sd(age_2020, na.rm = T),
            BMI_mean = mean(BMI, na.rm = T),
            BMI_sd = sd(BMI, na.rm = T),
            children_fathered_mean = mean(n_children_fathered, na.rm = T),
            children_fathered_sd = sd(n_children_fathered, na.rm = T),
            live_births_mean = mean(n_live_births, na.rm = T),
            live_births_sd = sd(n_live_births, na.rm = T),
            stillbirths_mean = mean(n_stillbirths, na.rm = T),
            stillbirths_sd = sd(n_stillbirths, na.rm = T),
            miscarr_mean = mean(n_miscarr, na.rm = T),
            miscarr_sd = sd(n_miscarr, na.rm = T),
            terminations_mean = mean(n_terminations, na.rm = T),
            terminations_sd = sd(n_terminations, na.rm = T),
            age_menarche_mean = mean(age_menarche, na.rm = T),
            age_menarche_sd = sd(age_menarche, na.rm = T),
            age_menopause_mean = mean(age_menopause, na.rm = T),
            age_menopause_sd = sd(age_menopause, na.rm = T))
quant_summary <- t(quant_summary)

d %>% group_by(sex, past_smoker) %>% summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

d %>% group_by(sex, current_smoker) %>% summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

d %>% group_by(sex, alc_intake) %>% summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

d %>% group_by(sex, SSMT) %>% summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

d %>% group_by(sex, ever_oral_contraceptive) %>% summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

d %>% group_by(sex, ever_HRT) %>% summarise(n = n()) %>% 
  mutate(freq = n / sum(n))

d %>% group_by(sex, hysterectomy) %>% summarise(n = n()) %>% 
  mutate(freq = n / sum(n))