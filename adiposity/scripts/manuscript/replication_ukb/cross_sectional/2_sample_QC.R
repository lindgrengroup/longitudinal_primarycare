# Author: Samvida S. Venkatesh
# Date: 01/04/21

library(tidyverse)

# Read data ----

infile_path <- "" # REDACTED
ukb_path <- "" # REDACTED
outfile_path <- "" # REDACTED

STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")

trait_dat <- lapply(STRATA, function (st) {
  res <- read.table(paste0(infile_path, "/traits_for_gwas/", st, "_cross_sectional.txt"),
                      sep = "\t", header = T, stringsAsFactors = F)
    return (res)
})
names(trait_dat) <- STRATA
NPCs <- 21

# QC file from UKBB
qc <- read.table(paste0(ukb_path, "/QC/ukb_sqc_v2.txt"), header = T, 
                 na.string = c("NA", "", "."), stringsAsFactors = F)

# fam file corresponding to the QC file provided by UKBB
fam <- read.table(paste0(ukb_path, "/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"), 
                  header = F)
# Add IDs to QC file
qc$eid <- fam[, 1]

# Phenotype file from UKBB
pheno <- read.table(paste0(ukb_path, "/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv"),
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "eid"

# Prepare data for genotyping QC ----

for_gen_QC <- lapply(STRATA, function (st_name) {
  # Get QC file info for ids
  df <- qc %>% select(c("eid", "Submitted.Gender", "Inferred.Gender",
                        "het.missing.outliers", "excess.relatives",
                        "in.Phasing.Input.chr1_22", 
                        "in.white.British.ancestry.subset",
                        "putative.sex.chromosome.aneuploidy",
                        "sample.qc.missing.rate",
                        "in.kinship.table",
                        "excluded.from.kinship.inference",
                        "genotyping.array",
                        paste0("PC", 1:NPCs))) %>%
    filter(eid %in% trait_dat[[st_name]]$IID)
  # Merge phenotype file info 
  # (f.22001.0.0: genotyped and recommended exclusion)
  # (f.54.0.0: UKB assessment centre)
  df <- merge(df, pheno[, c("eid", "f.22001.0.0", "f.54.0.0")], by = "eid")
  colnames(df)[which(colnames(df) == "f.54.0.0")] <- "UKB_assmt_centre"
  return (df)
})
names(for_gen_QC) <- STRATA

# Genotyping QC functions ----

## Withdrawn consent ----

remove_withdrawn_ids <- function (data, qc_log_file) {
  
  # Path to UKBB provided list of individuals that have withdrawn consent
  withdrawn <- read.table(paste0(ukb_path, "/QC/w11867_20220525.csv"), 
                          header = F)
  cleaned <- subset(data, !(data$eid %in% withdrawn$V1))
  
  # Report QC metrics
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Individuals that withdrew consent: ", 
            nrow(data) - nrow(cleaned), "\n",
            "REMAINING, ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)
}

remove_negative_ids <- function (data, qc_log_file) {
  
  data$eid <- as.numeric(data$eid)
  cleaned <- subset(data, data$eid > 0)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Individuals with negative IDs (withdrawn consent): ", 
            nrow(data) - nrow(cleaned), "\n",
            "REMAINING, ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)
}

## Sex ---- 

qc_sex_mismatch <- function(data, qc_log_file) {
  
  cleaned <- subset(data, 
                    !is.na(data$Submitted.Gender) & !is.na(data$Inferred.Gender) & 
                      data$Submitted.Gender == data$Inferred.Gender)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Sex mismatch ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep=""))
  sink()
  
  return(cleaned)	
}	

## Ancestry ----

keep_white_british_ancestry <- function (data, qc_log_file) {
  
  cleaned <- subset(data, data$in.white.British.ancestry.subset == 1)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, Not in white British ancestry subset: ",
            length(which(data$in.white.British.ancestry.subset != 1)), "\n",
            "REMAINING: ",
            nrow(cleaned), "\n", sep = ""))
  sink()
  
  return (cleaned)
  
}

## Genotyping ----

qc_het_miss <- function (data, qc_log_file) {
  
  cleaned <- subset(data, !is.na(data$het.missing.outliers) & 
                      data$het.missing.outliers != 1)	
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, poor heterozygosity or missingness: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)
  
}

qc_not_in_phasing  <- function(data, qc_log_file) {
  
  cleaned <- subset(data, !is.na(data$in.Phasing.Input.chr1_22) & 
                      data$in.Phasing.Input.chr1_22 != 0)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, not used in autosome phasing: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)	
  
}		

qc_sex_chr_aneupl  <- function(data, qc_log_file) {
  
  cleaned <- subset(data, !is.na(data$putative.sex.chromosome.aneuploid) & 
                      data$putative.sex.chromosome.aneuploid != 1)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, putative sex chr aneuploidy: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)	
  
}	

## Relatedness ----

qc_excess_related <- function(data, qc_log_file) {
  
  cleaned <- subset(data, !is.na(data$excess.relatives) &
                      data$excess.relatives != 1)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, excess relatives (>10 3rd degree relatives): ",
            nrow(data) - nrow(cleaned),
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)
  
}

qc_related <- function(data, qc_log_file) {
  
  # Pathway to UKBB list of related individuals
  related <- read.table(paste0(ukb_path, "/QC/ukb1186_rel_s488366.dat"),
                        header = T)
  
  # For each pair of related individuals
  # remove the samples with the highest missingness
  related <- related[related$Kinship > 0.0884 &
                       related$ID1 %in% data$eid & related$ID2 %in% data$eid, ]
  
  related$miss1 = data$sample.qc.missing.rate[match(related$ID1, data$eid)]
  related$miss2 = data$sample.qc.missing.rate[match(related$ID2, data$eid)]
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
  
  cleaned <- subset(data, !(data$eid %in% related$id_remove))
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Relatedness pairs with errors: ",
            length(which(related$id_remove == "error")), "\n",
            "**FILTER** Individuals excluded because of relatedness: ",
            nrow(data[data$eid %in% related$id_remove, ]), "\n",
            "REMAINING NOT RELATED: ", nrow(cleaned), "\n\n", sep = ""))
  sink()
  
  return(cleaned)
  
}

qc_kinship_table <- function(data, qc_log_file) {
  
  cleaned <- subset(data, !is.na(data$excluded.from.kinship.inference) &
                      data$excluded.from.kinship.inference == 0)
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** Excluded from kinship inference: ",
            nrow(data) - nrow(cleaned),
            " ; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)
  
}

## Other exclusions ----

ukb_recommended_excl <- function (data, qc_log_file) {
  
  # Field: f.22010.0.0, coding: 1 - recommended exclusion
  
  cleaned <- data
  
  remove <- which(cleaned$f.22010.0.0 == 1)
  if (length(remove) > 0) { cleaned <- cleaned[-remove, ] }
  
  sink(qc_log_file, append = T)
  cat(paste("**FILTER** EXCLUDED, recommended UKBIOBANK exclusion: ", 
            nrow(data) - nrow(cleaned), 
            "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
  sink()
  
  return(cleaned)	
  
}

# Perform genotyping QC ----

qcd_ids <- lapply(STRATA, function (st_name) {
  data <- for_gen_QC[[st_name]]
  qc_log_file <- paste0(outfile_path, "/GWAS/log_files/", st_name, ".txt")
  
  # Print sample characteristics before QC
  sink(qc_log_file, append = T)
  cat(paste("TOTAL SAMPLE SIZE: ", nrow(data), "\n", sep = ""))
  cat(paste("   GENOTYPED ", length(which(!is.na(data$f.22001.0.0))), "\n", sep = ""))
  sink()
  
  # Remove individuals that have not been genotyped
  cleaned <- subset(data, !is.na(data$f.22001.0.0))
  
  # Apply QC functions
  cleaned <- remove_withdrawn_ids(cleaned, qc_log_file)
  cleaned <- remove_negative_ids(cleaned, qc_log_file)
  cleaned <- qc_sex_mismatch(cleaned, qc_log_file)
  cleaned <- keep_white_british_ancestry(cleaned, qc_log_file)
  cleaned <- qc_het_miss(cleaned, qc_log_file)
  cleaned <- qc_excess_related(cleaned, qc_log_file)
  # cleaned <- qc_related(cleaned, qc_log_file)
  # cleaned <- qc_kinship_table(cleaned, qc_log_file)
  cleaned <- ukb_recommended_excl(cleaned, qc_log_file)
  
  return (cleaned)
})
names(qcd_ids) <- STRATA

# Create id files for GWAS ----

ids_for_gwas <- lapply(STRATA, function (st_name) {
  to_write <- qcd_ids[[st_name]] %>% 
    mutate(FID = eid, IID = eid) %>%
    select(FID, IID, genotyping.array, paste0("PC", 1:NPCs))
  
  write.table(to_write, 
              paste0(outfile_path, "/GWAS/sample_qc/", 
                     st_name, "_ids_passed_qc.txt"), 
              sep = "\t", row.names = F, quote = F)
  return (to_write)
})
names(ids_for_gwas) <- STRATA

# Subset trait files to IDs that pass QC ----

write_traits <- lapply(STRATA, function (st_name) {
  keep_traits <- trait_dat[[st_name]] %>% 
    filter(IID %in% ids_for_gwas[[st_name]]$IID)
  
  write.table(keep_traits, 
              paste0(outfile_path, "/GWAS/traits_for_gwas/qcd_",
                     st_name, "_cross_sectional.txt"),
              sep = "\t", row.names = F, quote = F)
  
  return (keep_traits)
})

