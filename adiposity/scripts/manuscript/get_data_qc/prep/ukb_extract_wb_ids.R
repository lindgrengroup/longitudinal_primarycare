# Author: Samvida S. Venkatesh
# Date: 26/08/21

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED

# QC file from UKBB
qc <- read.table(paste0(mainpath, "/QC/ukb_sqc_v2.txt"), 
                 header = T, 
                 na.string = c("NA", "", "."), stringsAsFactors = F)

# fam file corresponding to the QC file provided by UKBB
fam <- read.table(paste0(mainpath, "/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"), 
                  header = F)
# Add IDs to QC file
qc$eid <- fam[, 1]

# Remove withdrawn ids
withdrawn <- read.table(paste0(mainpath, "/QC/w11867_20210809.csv"), 
                        header = F)
qc <- subset(qc, !(qc$eid %in% withdrawn$V1))

wb_ids <- qc$eid[qc$in.white.British.ancestry.subset == 1]

write.table(wb_ids, 
            paste0(gen_resources_path, "/eids_white_british.txt"),
            sep = "\t", row.names = F, col.names = F, quote = F)
