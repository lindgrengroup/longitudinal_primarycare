# Author: Samvida S. Venkatesh
# Date: 01/12/2021

library(tidyverse)

# Read data ----

# Hormone "cases"
HORMONES <- read.table("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/hormone_list.txt")$V1
hormone_dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/data_popn_qcd_no_longit_filter.rds")[HORMONES]

# WB ancestry ids
WB_IDS <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/eids_white_british.txt",
                     sep = "\t", header = F, stringsAsFactors = F)$V1
# Annotated GP data
raw_gp_dat <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/gp_clinical_annotated.txt",
                         sep = "\t", header = T, comment.char = "$",
                         stringsAsFactors = F)

# Wrangle and manipulate data ----

# Retain WB subset of GP clinical
raw_gp_dat <- subset(raw_gp_dat, raw_gp_dat$eid %in% WB_IDS)

# Get IDs of all individuals that have had at least one hormone measured
# Remove these from control list
all_case_ids <- unique(unlist(lapply(hormone_dat, function (x) unique (x$eid))))
CONTROL_IDS <- WB_IDS[!WB_IDS %in% all_case_ids]
CONTROL_IDS <- data.frame(id = CONTROL_IDS,
                          sex = raw_gp_dat$sex[match(CONTROL_IDS,
                                                     raw_gp_dat$eid)],
                          dob = raw_gp_dat$dob[match(CONTROL_IDS,
                                                     raw_gp_dat$eid)])
CONTROL_IDS <- CONTROL_IDS[complete.cases(CONTROL_IDS), ]

write.table(CONTROL_IDS, 
            "/well/lindgren/UKBIOBANK/samvida/hormone_ehr/data/control_ids_no_hormones_measured.txt",
            sep = "\t", quote = F, row.names = F)
