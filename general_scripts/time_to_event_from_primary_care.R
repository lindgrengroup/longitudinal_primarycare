# Author: Samvida S. Venkatesh
# Date: 10/03/22

library(tidyverse)

# Read data ----

# Read dictionary to get to phenotype code-lists
dictionary <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F,
                         quote = "")

# Annotated GP data
gp_dat <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/gp_clinical_annotated.txt",
                          sep = "\t", header = T, comment.char = "$",
                          stringsAsFactors = F)

# Wrangle code data to only include codes in CPRD and exclude history ----

dictionary <- dictionary %>% filter(annot_CPRD != "")
UNIQ <- as.character(dictionary$unique_code)

code_lists <- lapply(UNIQ, function (ucode) {
  filename_to_get <- dictionary$annot_CPRD[dictionary$unique_code == ucode]
  codes <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/primary_care/",
                                          filename_to_get), 
                                   header = T, sep = "\t", stringsAsFactors = F,
                                   quote = "", comment.char = "~")
  codes <- subset(codes, !grepl("History", codes$Category))
  
  codes_V2 <- unique(codes$READV2_CODE)
  codes_V2 <- codes_V2[codes_V2 != "" & !is.na(codes_V2)]
  codes_V3 <- unique(codes$READV3_CODE)
  codes_V3 <- codes_V2[codes_V3 != "" & !is.na(codes_V3)]
  
  return (list(codes_V2 = codes_V2, codes_V3 = codes_V3))
})
names(code_lists) <- UNIQ

# Construct EID x phenotype matrix to get age at first recorded diagnosis ----

ALL_EIDS <- as.character(unique(gp_dat$eid))

sorted_gp_dat <- gp_dat %>% 
  group_by(eid) %>% 
  arrange(age_years, .by_group = T) %>%
  mutate(eid = as.character(eid)) %>%
  ungroup()

# Given disease V2 and V3 codes, function to return eids with age at event
# GP DATA HAS TO BE SORTED BY AGE BEFORE THIS FUNCTION CAN BE RUN
getFirstAge <- function (v2c, v3c) {
  df <- sorted_gp_dat %>% 
    filter(read_2 %in% v2c | read_3 %in% v3c) %>%
    group_by(eid) %>%
    summarise(age_at_first_diag = first(age_years))
  return (df)
}

eid_tte_lists <- lapply(UNIQ, function (ucode) {
  ids_with_tte <- getFirstAge(v2c = code_lists[[as.character(ucode)]]$codes_V2,
                              v3c = code_lists[[as.character(ucode)]]$codes_V3)
  colnames(ids_with_tte) <- c("eid", ucode)
  return (ids_with_tte)
})
eid_tte_matrix <- eid_tte_lists %>% reduce(full_join, by = "eid")

# Get time to event for first/last records in primary care
eid_tte_ends <- sorted_gp_dat %>% 
  group_by(eid) %>%
  summarise(age_at_first_record = first(age_years),
            age_at_last_record = last(age_years))

eid_tte_matrix <- full_join(eid_tte_ends, eid_tte_matrix, by = "eid")
colnames(eid_tte_matrix)[c(-1:-3)] <- UNIQ

# Save results
write.table(eid_tte_matrix, 
            "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)
