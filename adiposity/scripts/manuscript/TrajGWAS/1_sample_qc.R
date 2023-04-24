# Author: Samvida S. Venkatesh
# Date: 28/03/22

# NEED TO OUTPUT:
# Phenotype/covariates file
# File with row numbers of ids in bgen .sample file to include in analysis

library(tidyverse)

HORMONES <- c("FSH", "LH", "Oestradiol", "Progesterone", "Testosterone")

PCs <- paste0("PC", 1:21)
COVARS_FOR_MODEL <- c("sex", "data_provider", "age_event", 
                      "age_sq", "genotyping_array", PCs)

# Read data ----

# Hormone data (with multiple measures)
dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/gp_main_data_passed_longit_filter.rds")[HORMONES]

# General covariates
general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# QC file from UKBB
qc <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb_sqc_v2.txt", 
                 header = T, na.string = c("NA", "", "."), stringsAsFactors = F)
# fam file corresponding to the QC file provided by UKBB
fam <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam", 
                  header = F)
# Add IDs to QC file
qc$eid <- fam[, 1]

# Withdrawn
withdrawn <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/w11867_20220222.csv", 
                        header = F)$V1
withdrawn <- as.character(withdrawn)

# List of related individuals 
related <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb1186_rel_s488366.dat",
                      header = T)
related <- related %>% mutate(ID1 = as.character(ID1),
                              ID2 = as.character(ID2))

# Prepare data for QC ----

qc <- qc %>% select(c("eid", "Submitted.Gender", "Inferred.Gender",
                      "het.missing.outliers", "excess.relatives",
                      "in.Phasing.Input.chr1_22", 
                      "in.white.British.ancestry.subset",
                      "putative.sex.chromosome.aneuploidy",
                      "sample.qc.missing.rate",
                      "in.kinship.table",
                      "excluded.from.kinship.inference",
                      "genotyping.array",
                      PCs)) %>%
  mutate(eid = as.character(eid))

# QC functions ----
## Withdrawn consent ----

remove_withdrawn_ids <- function (ids, qc_log_file) {
  cleaned <- ids[!ids %in% withdrawn]
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Individuals that withdrew consent: ", 
            length(ids) - length(cleaned), "\n",
            "REMAINING, ", length(cleaned), "\n", sep = ""))
  sink()
  cleaned <- as.character(cleaned)
  return (cleaned)
}

remove_negative_ids <- function (ids, qc_log_file) {
  ids <- as.numeric(ids)
  cleaned <- ids[ids > 0]
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Individuals with negative IDs (withdrawn consent): ", 
            length(ids) - length(cleaned), "\n",
            "REMAINING, ", length(cleaned), "\n", sep = ""))
  sink()
  cleaned <- as.character(cleaned)
  return (cleaned)
}

## Sex ---- 

qc_sex_mismatch <- function(ids, qc_log_file) {
  data <- qc %>% filter(eid %in% ids)
  cleaned <- subset(data, 
                    !is.na(data$Submitted.Gender) & !is.na(data$Inferred.Gender) & 
                      data$Submitted.Gender == data$Inferred.Gender)
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Sex mismatch ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep=""))
  sink()
  return (cleaned$eid)	
}	

## Ancestry ----

keep_white_british_ancestry <- function (ids, qc_log_file) {
  data <- qc %>% filter(eid %in% ids)
  cleaned <- subset(data, data$in.white.British.ancestry.subset == 1)
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, Not in white British ancestry subset: ",
            length(which(data$in.white.British.ancestry.subset != 1)), "\n",
            "REMAINING: ",
            nrow(cleaned), "\n", sep = ""))
  sink()
  return (cleaned$eid)
}

## Genotyping ----

qc_het_miss <- function (ids, qc_log_file) {
  data <- qc %>% filter(eid %in% ids)
  cleaned <- subset(data, !is.na(data$het.missing.outliers) & 
                      data$het.missing.outliers != 1)	
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, poor heterozygosity or missingness: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  return(cleaned$eid)
  
}

qc_not_in_phasing  <- function(ids, qc_log_file) {
  data <- qc %>% filter(eid %in% ids)
  cleaned <- subset(data, !is.na(data$in.Phasing.Input.chr1_22) & 
                      data$in.Phasing.Input.chr1_22 != 0)
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, not used in autosome phasing: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  return(cleaned$eid)	
}		

qc_sex_chr_aneupl  <- function(ids, qc_log_file) {
  data <- qc %>% filter(eid %in% ids)
  cleaned <- subset(data, !is.na(data$putative.sex.chromosome.aneuploid) & 
                      data$putative.sex.chromosome.aneuploid != 1)
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, putative sex chr aneuploidy: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  return(cleaned$eid)	
}	

## Relatedness ----

qc_excess_related <- function(ids, qc_log_file) {
  data <- qc %>% filter(eid %in% ids)
  cleaned <- subset(data, !is.na(data$excess.relatives) &
                      data$excess.relatives != 1)
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, excess relatives (>10 3rd degree relatives): ",
            nrow(data) - nrow(cleaned),
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  return(cleaned$eid)
}

qc_related <- function(ids, qc_log_file) {
  
  data <- qc %>% filter(eid %in% ids)
  
  # For each pair of related individuals
  # remove the samples with the highest missingness
  related <- related[related$Kinship > 0.0884 &
                       related$ID1 %in% ids & related$ID2 %in% ids, ]
  
  related$miss1 = data$sample.qc.missing.rate[match(related$ID1, ids)]
  related$miss2 = data$sample.qc.missing.rate[match(related$ID2, ids)]
  related$max_miss <- pmax(related$miss1, related$miss2)
  
  # Remove according to rule above
  related$id_remove <- ifelse(is.na(related$miss1) & is.na(related$miss2),
                              related$ID2,
                              ifelse(is.na(related$miss1), related$ID1,
                                     ifelse(is.na(related$miss2), related$ID2,
                                            ifelse(related$miss1 ==
                                                     related$max_miss, related$ID1,
                                                   ifelse(related$miss2 ==
                                                            related$max_miss,
                                                          related$ID2, "error")))))
  
  cleaned <- ids[!(ids %in% related$id_remove)]
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Relatedness pairs with errors: ",
            length(which(related$id_remove == "error")), "\n",
            "**FILTER** Individuals excluded because of relatedness: ",
            length(which(ids %in% related$id_remove)), "\n",
            "REMAINING NOT RELATED: ", length(cleaned), "\n\n", sep = ""))
  sink()
  return(cleaned)
}

qc_kinship_table <- function (ids, qc_log_file) {
  data <- qc %>% filter(eid %in% ids)
  cleaned <- subset(data, !is.na(data$excluded.from.kinship.inference) &
                      data$excluded.from.kinship.inference == 0)
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Excluded from kinship inference: ",
            nrow(data) - nrow(cleaned),
            " ; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  return (cleaned$eid)
}

# Perform genotyping QC ----

qcd_ids <- lapply(HORMONES, function (hr) {
  ids_to_qc <- unique(dat[[hr]]$eid)
  
  qc_log_file <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/TrajGWAS/log_files/sample_qc_",
                        hr, ".txt")
  
  # Print sample characteristics before QC
  sink(qc_log_file, append = T)
  cat(paste("TOTAL SAMPLE SIZE: ", length(ids_to_qc), "\n", sep = ""))
  sink()
  
  # Apply QC functions
  cleaned <- remove_withdrawn_ids(ids_to_qc, qc_log_file)
  cleaned <- remove_negative_ids(cleaned, qc_log_file)
  cleaned <- qc_sex_mismatch(cleaned, qc_log_file)
  cleaned <- keep_white_british_ancestry(cleaned, qc_log_file)
  cleaned <- qc_het_miss(cleaned, qc_log_file)
  cleaned <- qc_excess_related(cleaned, qc_log_file)
  cleaned <- qc_related(cleaned, qc_log_file)
  cleaned <- qc_kinship_table(cleaned, qc_log_file)
  
  return (cleaned)
})
names(qcd_ids) <- HORMONES

# Write covariates + phenotypes file for each sex ----

write_dat <- lapply(HORMONES, function (hr) {
  pheno_keep <- dat[[hr]] %>% 
    filter(eid %in% qcd_ids[[hr]]) %>%
    select(all_of(c("eid", "value", "data_provider", "age_event"))) %>%
    mutate(eid = as.character(eid),
           age_sq = age_event^2)
  
  res <- left_join(pheno_keep, 
                           general_covars[, c("eid", "sex", PCs)],
                           by = "eid") %>%
    mutate(genotyping_array = qc$genotyping.array[match(eid, qc$eid)])
  
  res <- res[, c("eid", "value", COVARS_FOR_MODEL)]
  
  # Write male, female, sex-combined bgen sample row ids
  sex_strata <- c("F", "M", "sex_comb")
  for (ss in sex_strata) {
    to_write <- res
    if (ss != "sex_comb") to_write <- to_write %>% filter(sex == ss)
    
    write.csv(to_write, 
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/TrajGWAS/data/",
                     hr, "_", ss, "_pheno_covars.csv"),
              row.names = F, quote = F)
  }
  return ()
})

