# Author: Samvida S. Venkatesh
# Date: 03/05/22

library(survival)
library(survminer)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
theme_set(theme_bw())

# Read and wrangle data ----

survival_log <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/qc/survival_data_removed_individuals.txt"

# Full UKB phenotypes file to get sex
# Main UKBB phenotypes file
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

annot_dictionary <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                               sep = "\t", header = T, stringsAsFactors = F, 
                               quote = "", comment.char = "$")
annot_dictionary$phenotype <- gsub('"', "", annot_dictionary$phenotype)

DIAGNOSES <- annot_dictionary$phenotype

age_at_diag_matrix <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt",
                                 sep = "\t", header = T, stringsAsFactors = F,
                                 quote = "")
colnames(age_at_diag_matrix) <- gsub("^X", "", colnames(age_at_diag_matrix))
colnames(age_at_diag_matrix)[c(-1:-5)] <- DIAGNOSES

age_at_diag_matrix$death <- general_covars$age_at_death[match(age_at_diag_matrix$eid,
                                                              general_covars$eid)]

age_at_diag_matrix <- age_at_diag_matrix %>%
  mutate(across(all_of(c("age_at_first_GP_record",
                         "age_at_last_GP_record",
                         "age_at_first_HES_record", "age_at_last_HES_record",
                         "death")), 
                as.numeric)) %>%
  mutate(eid = as.character(eid))

# Create survival information ----

# Get age at which observations first started (start of window)
# and age at which there are no longer data (end of window)
# If death recorded after last HES/GP event, include that in the censoring info
# assuming that anything else between the last HES/GP event and death is 
# not relevant to the diagnosis
# Make sure we have an observation "window" and not just one event
age_at_diag_matrix <- age_at_diag_matrix %>%
  mutate(age_at_first_record = pmin(age_at_first_GP_record,
                                    age_at_first_HES_record, 
                                    na.rm = T),
         age_at_last_record = pmax(age_at_last_GP_record,
                                   age_at_last_HES_record,
                                   na.rm = T)) %>%
  filter(!is.na(age_at_first_record) & !is.na(age_at_last_record)
         & age_at_first_record < age_at_last_record)

sink(survival_log, append = T)
cat(paste0("Number of individuals with non-conflicting survival information: ", 
           nrow(age_at_diag_matrix), "\n", "\n"))
sink()

# For each diagnosis, create a df with survival time and censoring info
surv_dat <- lapply(c("death", DIAGNOSES), function (diag) {
  
  res <- age_at_diag_matrix %>% 
    select(all_of(c("eid", "age_at_first_record", "age_at_last_record",
                    diag))) %>%
    rename(age_at_event = !!as.symbol(diag))
  
  sink(survival_log, append = T)
  cat(paste0("Diagnosis: ", diag, "\n",
             "\t", "Number of cases: ", sum(!is.na(res$age_at_event)), "\n",
             "\t", "Cases excluded with history prior to observation window: ", 
             sum(res$age_at_event == "history", na.rm = T), "\n"))
  sink()
  
  # Remove individuals with history prior to observation
  have_hist <- which(res$age_at_event == "history")
  cleaned <- res
  if (length(have_hist) > 0) {
    cleaned <- res[-have_hist, ]
  } 
  # If non-NA age at event, then event is observed - if not, then censored
  cleaned <- cleaned %>% 
    mutate(event_observed = !is.na(age_at_event)) 
  cleaned$age_at_event[!cleaned$event_observed] <- 
    cleaned$age_at_last_record[!cleaned$event_observed]
  cleaned$age_at_event <- as.numeric(cleaned$age_at_event)
  
  # Remove any inconsistencies, i.e. when age at disease onset is 
  # earlier than the first record in GP or HES
  # If age at disease is itself the first record in GP or HES,
  # tweak the first record to reduce age by 0.01 yrs for 
  # computational ease
  tweak_age <- which(cleaned$age_at_event == cleaned$age_at_first_record)
  cleaned$age_at_first_record[tweak_age] <- 
    cleaned$age_at_first_record[tweak_age] - 0.01
  
  cleaned <- cleaned %>% filter(age_at_event > age_at_first_record)
  sink(survival_log, append = T)
  cat(paste0("\t", "FINAL cases with non-conflicting survival data: ",  
             sum(cleaned$event_observed), "\n"))
  sink()
  
  return (cleaned)
})
names(surv_dat) <- c("death", DIAGNOSES)

# Sanity check by plotting survival curves for full population (sex)
# COMMENT THIS OUT ONCE DONE

# Add covariates
surv_dat_sex <- lapply(surv_dat, function (df) {
  df$eid <- as.character(df$eid)
  df$sex <- general_covars$sex[match(df$eid, general_covars$eid)]
  df <- df %>% filter(!is.na(sex))
  return (df)
})

# Plot Kaplan Meier curves stratified by sex
km_plots <- lapply(c("death", DIAGNOSES), function (diag) {
  df <- surv_dat_sex[[diag]]
  fitted_dat <- survfit(Surv(age_at_first_record,
                             age_at_event,
                             event_observed) ~ sex, data = df)
  res_plot <- ggsurvplot_df(surv_summary(fitted_dat),
                            fun = "event",
                            conf.int = T, censor = F,
                            ggtheme = theme_bw(),
                            xlab = "Age (years)", ylab = "Cumulative incidence (%)",
                            title = diag)
  diag_name_tmp <- gsub("/", "_", diag)
  png(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/plots/kaplan_meier_survival/",
             diag_name_tmp, ".png"))
  print(res_plot)
  dev.off()
  return ()
})

# This looks reasonable so we can save survival data for future use
saveRDS(surv_dat,
        "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_survival_dat.rds")
