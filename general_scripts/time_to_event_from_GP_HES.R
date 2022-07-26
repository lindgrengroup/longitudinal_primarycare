# Author: Samvida S. Venkatesh
# Date: 10/03/22

library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

# Read data ----

# Read phenotype code-lists
dictionary <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F,
                         na.string = c("NA", "", "."), quote = "")
dictionary$phenotype <- gsub('\\"', "", dictionary$phenotype)
dictionary$ICD_chapter_desc <- gsub('\\"', "", dictionary$ICD_chapter_desc)

# Annotated GP data
gp_dat <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/gp_clinical_annotated.txt",
                     sep = "\t", header = T, comment.char = "$",
                     stringsAsFactors = F)
gp_dat <- gp_dat %>% 
  mutate(across(all_of(c("eid", "read_2", "read_3")), as.character))

# Annotated HES data
raw_hes_dat <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/hes_age_annotated.txt",
                          sep = "\t", header = T, comment.char = "~", 
                          stringsAsFactors = F, quote = "",
                          na.string = c("NA", "", "."))
# Subset to ICD- and OPCS-relevant columns
ICD9COLS <- colnames(raw_hes_dat)[grep("icd9", colnames(raw_hes_dat))]
ICD10COLS <- colnames(raw_hes_dat)[grep("icd10", colnames(raw_hes_dat))]
OPCSCOLS <- colnames(raw_hes_dat)[grep("oper4", colnames(raw_hes_dat))]

hes_dat <- raw_hes_dat %>% 
  select(all_of(c("eid", "record_id", "age_record", 
                  ICD9COLS, ICD10COLS, OPCSCOLS))) %>%
  mutate(across(all_of(c("eid", ICD9COLS, ICD10COLS, OPCSCOLS)), 
                as.character)) 

# Annotated death data
raw_death_dat <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/death_register_age_annotated.txt",
                            sep = "\t", header = T, comment.char = "~", 
                            stringsAsFactors = F, quote = "",
                            na.string = c("NA", "", "."))
death_dat <- raw_death_dat %>%
  select(all_of(c("eid", "age_at_death", "cause_icd10")))  %>%
  mutate(across(all_of(c("eid", "cause_icd10")), 
                as.character)) 

# Wrangle code data to separate diagnosis and history codes ----

UNIQ <- as.character(dictionary$unique_code)

code_lists <- lapply(UNIQ, function (ucode) {
  # Get CPRD codes
  cprd_file <- dictionary$annot_CPRD[dictionary$unique_code == ucode]
  if (!is.na(cprd_file)) {
    codes <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/primary_care/",
                               cprd_file), 
                        header = T, sep = "\t", stringsAsFactors = F,
                        quote = "", comment.char = "~",
                        na.string = c("NA", "", "."))
    codes <- subset(codes, !grepl("History", codes$Category))
    
    codes_V2 <- unique(as.character(codes$READV2_CODE))
    codes_V2 <- codes_V2[!is.na(codes_V2)]
    codes_V3 <- unique(as.character(codes$READV3_CODE))
    codes_V3 <- codes_V2[!is.na(codes_V3)]
  } else {
    codes_V2 <- "NOMATCH"
    codes_V3 <- "NOMATCH"
  }
  
  # Get ICD9 and ICD10 codes
  icd_file <- dictionary$annot_ICD[dictionary$unique_code == ucode]
  if (!is.na(icd_file)) {
    codes <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/secondary_care/",
                               icd_file), 
                        header = T, sep = "\t", stringsAsFactors = F,
                        quote = "", comment.char = "~",
                        na.string = c("NA", "", "."))
    codes <- subset(codes, !grepl("History", codes$Category))
    
    codes_ICD10 <- unique(as.character(codes$ICD10))
    codes_ICD10 <- codes_ICD10[!is.na(codes_ICD10) & codes_ICD10 != "UNDEF"]
    codes_ICD9 <- unique(as.character(codes$ICD9))
    codes_ICD9 <- codes_ICD9[!is.na(codes_ICD9) & codes_ICD9 != "UNDEF"]
  } else {
    codes_ICD10 <- "NOMATCH"
    codes_ICD9 <- "NOMATCH"
  }
  
  # Get OPCS4 codes
  opcs_file <- dictionary$OPCS[dictionary$unique_code == ucode]
  if (!is.na(opcs_file)) {
    codes <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/secondary_care/",
                               opcs_file), 
                        header = T, sep = ",", stringsAsFactors = F,
                        comment.char = "~",
                        na.string = c("NA", "", "."))
    codes <- subset(codes, !grepl("History", codes$Category))
    
    codes_OPCS <- gsub("\\.", "", codes$OPCS4code)
    codes_OPCS <- unique(as.character(codes_OPCS))
    codes_OPCS <- codes_OPCS[!is.na(codes_OPCS)]
  } else {
    codes_OPCS <- "NOMATCH"
  }
  
  return (list(codes_V2 = codes_V2, codes_V3 = codes_V3,
               codes_ICD9 = codes_ICD9, codes_ICD10 = codes_ICD10,
               codes_OPCS = codes_OPCS))
})
names(code_lists) <- UNIQ

history_code_lists <- lapply(UNIQ, function (ucode) {
  # Get CPRD codes
  cprd_file <- dictionary$annot_CPRD[dictionary$unique_code == ucode]
  if (!is.na(cprd_file)) {
    codes <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/primary_care/",
                               cprd_file), 
                        header = T, sep = "\t", stringsAsFactors = F,
                        quote = "", comment.char = "~",
                        na.string = c("NA", "", "."))
    codes <- subset(codes, grepl("History", codes$Category))
    
    codes_V2 <- unique(as.character(codes$READV2_CODE))
    codes_V2 <- codes_V2[!is.na(codes_V2)]
    codes_V3 <- unique(as.character(codes$READV3_CODE))
    codes_V3 <- codes_V2[!is.na(codes_V3)]
  } else {
    codes_V2 <- "NOMATCH"
    codes_V3 <- "NOMATCH"
  }
  
  # Get ICD9 and ICD10 codes
  icd_file <- dictionary$annot_ICD[dictionary$unique_code == ucode]
  if (!is.na(icd_file)) {
    codes <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/secondary_care/",
                               icd_file), 
                        header = T, sep = "\t", stringsAsFactors = F,
                        quote = "", comment.char = "~",
                        na.string = c("NA", "", "."))
    codes <- subset(codes, grepl("History", codes$Category))
    
    codes_ICD10 <- unique(as.character(codes$ICD10))
    codes_ICD10 <- codes_ICD10[!is.na(codes_ICD10) & codes_ICD10 != "UNDEF"]
    codes_ICD9 <- unique(as.character(codes$ICD9))
    codes_ICD9 <- codes_ICD9[!is.na(codes_ICD9) & codes_ICD9 != "UNDEF"]
  } else {
    codes_ICD10 <- "NOMATCH"
    codes_ICD9 <- "NOMATCH"
  }
  
  # Get OPCS4 codes
  opcs_file <- dictionary$OPCS[dictionary$unique_code == ucode]
  if (!is.na(opcs_file)) {
    codes <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/secondary_care/",
                               opcs_file), 
                        header = T, sep = ",", stringsAsFactors = F,
                        comment.char = "~",
                        na.string = c("NA", "", "."))
    codes <- subset(codes, grepl("History", codes$Category))
    
    codes_OPCS <- gsub("\\.", "", codes$OPCS4code)
    codes_OPCS <- unique(as.character(codes_OPCS))
    codes_OPCS <- codes_OPCS[!is.na(codes_OPCS)]
  } else {
    codes_OPCS <- "NOMATCH"
  }
  
  return (list(hist_codes_V2 = codes_V2, hist_codes_V3 = codes_V3,
               hist_codes_ICD9 = codes_ICD9, hist_codes_ICD10 = codes_ICD10,
               hist_codes_OPCS = codes_OPCS))
})
names(history_code_lists) <- UNIQ

# Wrangle GP, HES, and death data to arrange by event dates ----

sorted_gp_dat <- gp_dat %>% 
  group_by(eid) %>% 
  arrange(age_years, .by_group = T) %>% 
  ungroup()

sorted_hes_dat <- hes_dat %>% 
  group_by(eid) %>%
  arrange(age_record, .by_group = T) %>%
  ungroup()

# need to do this because some people have multiple death records
# with multiple causes of death recorded
sorted_death_dat <- death_dat %>%
  group_by(eid) %>%
  arrange(age_at_death, .by_group = T) %>%
  ungroup()

# Functions to get age at first diagnosis for specific codes ----

# Given disease V2 and V3 codes, function to return eids with age at event
# GP DATA HAS TO BE SORTED BY AGE BEFORE THIS FUNCTION CAN BE RUN
getFirstAgeGP <- function (v2c, v3c) {
  df <- sorted_gp_dat %>% 
    filter(read_2 %in% v2c | read_3 %in% v3c) %>%
    group_by(eid) %>%
    summarise(age_at_first_GP = first(age_years))
  return (df)
}

# Given disease ICD9, ICD10, or OPCS codes, 
# function to return eids with age at event
# HES DATA HAS TO BE SORTED BY AGE BEFORE THIS FUNCTION CAN BE RUN

getFirstAgeHES <- function (icd9c, icd10c, opcsc) {
  # First get case status
  df <- sorted_hes_dat %>% 
    mutate(across(all_of(ICD9COLS), ~ .x %in% icd9c, .names = "case_{.col}"),
           across(all_of(ICD10COLS), ~ .x %in% icd10c, .names = "case_{.col}"),
           across(all_of(OPCSCOLS), ~ .x %in% opcsc, .names = "case_{.col}"))
  df$case_status <- rowSums(df[, grep("^case_", colnames(df))], na.rm = T)
  
  # Then get first time the case was recorded
  res <- df %>% filter(case_status > 0) %>%
    group_by(eid) %>%
    summarise(age_at_first_HES = first(age_record))
  
  return (res)
}

# Given disease ICD10 codes,
# function to return eids with age at death caused by this event
# DEATH DATA HAS TO BE SORTED BY AGE BEFORE THIS FUNCTION CAN BE RUN
getFirstAgeDeath <- function (icd10c) {
  df <- sorted_death_dat %>% 
    filter(cause_icd10 %in% icd10c) %>%
    group_by(eid) %>%
    summarise(age_at_death = first(age_at_death))
  return (df)
}

# Given ages at history codes or diagnosis codes,
# ensure these align by marking "history" for any sets of ages where
# history precedes diagnosis, or where there is only a history and no diagnosis
combineDiagnosisHistoryInfo <- function (hist_age, diag_age) {
  if (is.na(hist_age) & is.na(diag_age)) res <- NA
  # only history present
  else if (!is.na(hist_age) & is.na(diag_age)) res <- "history"
  # both history and diagnosis present, but history before diagnosis
  # which makes it inconsistent
  else if (!is.na(hist_age) & !is.na(diag_age) & hist_age < diag_age) 
    res <- "history"
  else
    res <- diag_age
  return (res)
}

# Apply functions and combine GP, HES, death register results ----

eid_tte_lists <- lapply(UNIQ, function (ucode) {
  ids_with_gp_record <- 
    getFirstAgeGP(v2c = code_lists[[as.character(ucode)]]$codes_V2,
                  v3c = code_lists[[as.character(ucode)]]$codes_V3)
  
  ids_with_hes_record <- 
    getFirstAgeHES(icd9c = code_lists[[as.character(ucode)]]$codes_ICD9,
                   icd10c = code_lists[[as.character(ucode)]]$codes_ICD10,
                   opcsc = code_lists[[as.character(ucode)]]$codes_OPCS)
  
  ids_with_death_record <- 
    getFirstAgeDeath(icd10c = code_lists[[as.character(ucode)]]$codes_ICD10)
  
  ids_with_tte <- full_join(ids_with_gp_record, ids_with_hes_record, 
                            by = "eid") %>%
    full_join(ids_with_death_record, by = "eid") %>%
    mutate(age_at_first_diag = pmin(age_at_first_GP, age_at_first_HES, age_at_death,
                                    na.rm = T)) %>%
    select(all_of(c("eid", "age_at_first_diag")))
  
  # Get individuals who have history of a diagnosis and the age at which 
  # history was recorded
  # This is not relevant for the death codes
  ids_with_gp_history <- 
    getFirstAgeGP(v2c = history_code_lists[[as.character(ucode)]]$codes_V2,
                  v3c = history_code_lists[[as.character(ucode)]]$codes_V3)
  
  ids_with_hes_history <- 
    getFirstAgeHES(icd9c = history_code_lists[[as.character(ucode)]]$codes_ICD9,
                   icd10c = history_code_lists[[as.character(ucode)]]$codes_ICD10,
                   opcsc = history_code_lists[[as.character(ucode)]]$codes_OPCS)
  
  ids_with_history <- full_join(ids_with_gp_history, ids_with_hes_history, 
                                by = "eid") %>%
  mutate(age_at_first_history = pmin(age_at_first_GP, age_at_first_HES, 
                                     na.rm = T)) %>%
    select(all_of(c("eid", "age_at_first_history")))
  
  # Combine the diagnosis and history data
  # If there is a history code preceding a diagnosis code, then the individual
  # has to be censored from survival data, so change age-at-diagnosis to "history"
  full_dat <- full_join(ids_with_tte, ids_with_history,
                        by = "eid") %>%
    rowwise() %>%
    mutate(age_at_first_diag = combineDiagnosisHistoryInfo(age_at_first_history,
                                                           age_at_first_diag))

  # Return age at diagnosis for matrix
  for_mat <- full_dat %>% 
    select(all_of(c("eid", "age_at_first_diag")))
  colnames(for_mat) <- c("eid", ucode)
  
  return (for_mat)
})
names(eid_tte_lists) <- UNIQ

# Build eid-time-to-event matrix ----

eid_tte_matrix <- eid_tte_lists %>% reduce(full_join, by = "eid")

# Get time to event for first/last records in primary care
eid_tte_GP_ends <- sorted_gp_dat %>% 
  group_by(eid) %>%
  summarise(age_at_first_GP_record = first(age_years),
            age_at_last_GP_record = last(age_years))

eid_tte_HES_ends <- sorted_hes_dat %>% 
  group_by(eid) %>%
  summarise(age_at_first_HES_record = first(age_record),
            age_at_last_HES_record = last(age_record))

# Get age at death (first)
eid_age_death <- sorted_death_dat %>%
  group_by(eid) %>%
  summarise(age_at_death = first(age_at_death))

eid_tte_ends <- full_join(eid_tte_GP_ends, 
                          eid_tte_HES_ends, by = "eid") %>%
  full_join(eid_age_death, by = "eid")

eid_tte_matrix <- full_join(eid_tte_ends, 
                            eid_tte_matrix, by = "eid")
colnames(eid_tte_matrix)[c(-1:-6)] <- UNIQ

# Save results
write.table(eid_tte_matrix, 
            "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)

# # Plot summary statistics on GP vs HES ----
# 
# # Get data to plot
# prop_cases_gp_hes <- lapply(UNIQ, function (ucode) {
#   return (eid_tte_lists[[as.character(ucode)]]$prop_summary)
# })
# prop_cases_gp_hes <- bind_rows(prop_cases_gp_hes)
# saveRDS(prop_cases_gp_hes,
#         "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_tte_prop_gp_hes_both.rds")
# 
# prop_first_gp <- lapply(UNIQ, function (ucode) {
#   return (eid_tte_lists[[as.character(ucode)]]$first_summary)
# })
# prop_first_gp <- bind_rows(prop_first_gp)
# saveRDS(prop_first_gp,
#         "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_tte_prop_first_gp.rds")
# 
# ## Generic plotting tools ----
# 
# ALL_CHAPS <- sort(unique(dictionary$ICD_chapter))
# DIAG_LEVELS <- dictionary$phenotype[order(dictionary$ICD_chapter)]
# ICD_DESC_LEVELS <- unique(dictionary$ICD_chapter_desc[order(dictionary$ICD_chapter)])
# 
# REC_LEVELS <- c("primary_care", "secondary_care", "both")
# col_palette <- c("#E41A1C", "#377EB8", "#4DAF4A")
# names(col_palette) <- REC_LEVELS
# 
# # Mutate df as necessary
# wrangleWithDict <- function (df, var_to_plot) {
#   diag_res <- df %>%
#     mutate(rec_from = factor(!!as.symbol(var_to_plot), levels = REC_LEVELS),
#            diagnosis = factor(dictionary$phenotype[match(ucode,
#                                                          dictionary$unique_code)],
#                               levels = DIAG_LEVELS),
#            ICD_chapter = factor(dictionary$ICD_chapter[match(diagnosis,
#                                                              dictionary$phenotype)],
#                                 levels = ALL_CHAPS),
#            ICD_chapter_desc = factor(dictionary$ICD_chapter_desc[match(diagnosis,
#                                                                        dictionary$phenotype)],
#                                      levels = ICD_DESC_LEVELS))
#   
#   # Summarise by ICD chapter
#   chap_res <- diag_res %>%
#     group_by(ICD_chapter_desc) %>% mutate(n_chapter = sum(n)) %>%
#     ungroup() %>%
#     group_by(ICD_chapter_desc, !!as.symbol(var_to_plot)) %>%
#     summarise(n_chap_per_var = sum(n),
#               mean_prop = n_chap_per_var / n_chapter) %>%
#     distinct()
#   
#   return (list(diagnosis_res = diag_res,
#                chapter_res = chap_res))
# }
# 
# # Return diagnosis-level and chapter-level plots
# getPlots <- function (df_list, 
#                       var_to_plot, legend_label) {
#   
#   # Split this by disease chapter into multiple plots
#   per_disease <- lapply(ALL_CHAPS, function (cfilter) {
#     to_plot <- df_list$diagnosis_res %>% filter(ICD_chapter == cfilter)
#     if (nrow(to_plot) == 0) {
#       res <- NULL
#     } else {
#       res <-  ggplot(to_plot, aes(x = diagnosis, y = n, 
#                                   fill = !!as.symbol(var_to_plot))) +
#         geom_bar(position = "stack", stat = "identity") + 
#         scale_x_discrete(label = function(x) substr(x, 1, 20)) +
#         scale_fill_manual(values = col_palette) +
#         labs(x = "diagnosis", y = "number of cases", fill = legend_label) +
#         coord_flip()
#     }
#     return (res)
#   })
#   # Combine into multiple plots on multiple pages
#   per_disease_return_print <- ggarrange(plotlist = per_disease, 
#                                         nrow = 2, ncol = 1,
#                                         common.legend = T)
#   
#   per_chap <- ggplot(df_list$chapter_res, 
#                      aes(x = ICD_chapter_desc, y = mean_prop,
#                          fill = !!as.symbol(var_to_plot))) +
#     geom_bar(position = "stack", stat = "identity") + 
#     scale_x_discrete(label = function(x) substr(x, 1, 20)) +
#     scale_fill_manual(values = col_palette) +
#     labs(x = "ICD_chapter", y = "proportion of cases", fill = legend_label) +
#     coord_flip()
#   
#   return (list(per_disease_return_print, per_chap))
# }
# 
# ## Apply plotting functions ----
# 
# # GP vs HES cases diagnosed
# diag_info <- wrangleWithDict(prop_cases_gp_hes, 
#                              var_to_plot = "rec_from")
# # GP vs HES which came first
# first_info <- wrangleWithDict(prop_first_gp,
#                               var_to_plot = "which_first")
# 
# pdf("plots/eid_tte_gp_vs_hes.pdf",
#     onefile = T)
# getPlots(diag_info, var_to_plot = "rec_from", legend_label = "recorded_in")
# getPlots(first_info, var_to_plot = "which_first", legend_label = "first_recorded_in")
# dev.off()
