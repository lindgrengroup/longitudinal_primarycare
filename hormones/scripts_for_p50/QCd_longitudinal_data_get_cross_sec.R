# Author: Samvida S. Venkatesh
# Date: 16/11/2021

library(tidyverse)

# Read data ----

dat <- readRDS("indiv_qcd_data.rds")

# QC log file
qc_log <- "cross_sectional_qc.txt"

# Get cross-sectional value for each individual ----

# Define cross-sectional value as the first observed value closest to
# an individual's median and get age at event for this observed value

cross_sec_dat <- dat %>% group_by(eid) %>%
  # Get median value and distance of each observation to median
  mutate(medn_value = median(value),
         dist_to_medn = abs(value - median(value))) %>%
  # Arrange by age to get first observation that is closest to median
  group_by(eid) %>% arrange(age_event, .by_group = T) %>%
  slice(which.min(dist_to_medn))  

saveRDS(cross_sec_dat, 
        "qcd_cross_sec_data.rds")