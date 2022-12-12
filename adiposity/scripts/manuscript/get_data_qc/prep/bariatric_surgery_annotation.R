# Author: Samvida S. Venkatesh
# Date: 16/02/21

library(tidyverse)
library(lubridate)

mainpath <- "" # REDACTED

# Read full primary care data
gp_clinical <- read.table(paste0(mainpath, "/gp_clinical_annotated.txt"),
                          sep = "\t", header = T, comment.char = "$")
gp_clinical$event_dt <- as.Date(gp_clinical$event_dt, "%Y-%m-%d")
gp_clinical$dob <- as.Date(gp_clinical$dob, "%Y-%m-%d")

# Read bariatric surgery codes
bariatric_codes <- read.table(paste0(mainpath, "/bariatric_surgery_codes.txt"),
                              sep = "\t", header = T, comment.char = "$")

BARCURRENT_CODES <- bariatric_codes$read_code[grepl("^Bariatric", 
                                                    bariatric_codes$category)]
BARHIST_CODES <- bariatric_codes$read_code[grepl("^History", 
                                                 bariatric_codes$category)]

bar_current_EIDDATES <- gp_clinical[gp_clinical$read_2 %in% 
                                      BARCURRENT_CODES |
                                      gp_clinical$read_3 %in% 
                                      BARCURRENT_CODES, c("eid", "event_dt")]
# Get earliest recorded surgery date
bar_current_EIDDATES <- bar_current_EIDDATES %>% group_by(eid) %>%
  summarise(event_dt = min(event_dt))

# Add individuals with a history of bariatric surgery for whom we don't 
# have dates of surgery

bar_hist_EIDS <- unique(gp_clinical$eid[gp_clinical$read_2 %in% 
                                          BARHIST_CODES |
                                          gp_clinical$read_3 %in% 
                                          BARHIST_CODES])
# Remove any IDs that are also recorded in "current" as this means we know
# the surgery date and can retain some records
bar_hist_EIDS <- bar_hist_EIDS[!bar_hist_EIDS %in% bar_current_EIDDATES$eid]

bar_EIDDATES <- bind_rows(bar_current_EIDDATES,
                                  data.frame(eid = bar_hist_EIDS, 
                                             event_dt = NA))

write.table(bar_EIDDATES, 
            paste0(mainpath, "/bariatric_surgery_records.txt"),
            sep = "\t", row.names = F, quote = F)