# Author: Samvida S. Venkatesh
# Date: 01/04/21
# ADAPTED FROM ADIPOSITY CHANGE CONSORTIUM SOP

library(tidyverse)

# Read data ----

# Slopes
final_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_adj_slopes.rds")
PHENOTYPES <- names(final_slopes)
STRATA <- names(final_slopes[[1]])

# QC file from UKBB
qc <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt", header = T, 
                 na.string = c("NA", "", "."), stringsAsFactors = F)

# fam file corresponding to the QC file provided by UKBB
fam <- read.table("/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam", 
                  header = F)
# Add IDs to QC file
qc$eid <- fam[, 1]

# Phenotype file from UKBB
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "eid"

# Prepare data for genotyping QC ----

GWAS_STRATA <- STRATA[grep("white", STRATA)]

for_gen_QC <- lapply(PHENOTYPES, function (p) {
  res <- lapply(GWAS_STRATA, function (s) {
    df <- final_slopes[[p]][[s]]
    # Merge QC file info with slopes
    df <- merge(df, qc[, c("eid", "Submitted.Gender", "Inferred.Gender",
                           "het.missing.outliers", "excess.relatives",
                           "in.Phasing.Input.chr1_22", 
                           "in.white.British.ancestry.subset",
                           "putative.sex.chromosome.aneuploidy",
                           "sample.qc.missing.rate",
                           "in.kinship.table",
                           "excluded.from.kinship.inference",
                           "genotyping.array")], by = "eid")
    # Merge phenotype file info (f.22001.0.0: genotyped and recommended exclusion)
    df <- merge(df, pheno[, c("eid", "f.22001.0.0")], by = "eid")
    return (df)
  })
  names(res) <- GWAS_STRATA
  return (res)
})
names(for_gen_QC) <- PHENOTYPES

# Genotyping QC functions ----

## Withdrawn consent ----

remove_withdrawn_ids <- function (data, qc_log_file) {
  
  # Path to UKBB provided list of individuals that have withdrawn consent
  withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20210201.csv", 
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

# ## Relatedness ----
# 
# # DON'T HAVE TO PERFORM RELATEDNESS QC AS SUGEN ACCOUNTS FOR THIS
# # MAY HAVE TO REMOVE FOR OTHER METHODS
# 
# qc_excess_related <- function(data, qc_log_file) {
#   
#   cleaned <- subset(data, !is.na(data$excess.relatives) & 
#                       data$excess.relatives != 1)
#   
#   sink(qc_log_file, append = T)
#   cat(paste("**FILTER** EXCLUDED, excess relatives (>10 3rd degree relatives): ", 
#             nrow(data) - nrow(cleaned), 
#             "; REMAINING: ", nrow(cleaned), "\n", sep = ""))
#   sink()
#   
#   return(cleaned)	
#   
# }	
# 
# qc_related <- function(data, qc_log_file) {
# 
#   # Pathway to UKBB list of related individuals
#   related <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb1186_rel_s488366.dat",
#                         header = T)
# 
#   # For each pair of related individuals
#   # remove the samples with the highest missingness
#   related <- related[related$Kinship > 0.0884 &
#                        related$ID1 %in% data$eid & related$ID2 %in% data$eid, ]
# 
#   related$miss1 = data$sample.qc.missing.rate[match(related$ID1, data$eid)]
#   related$miss2 = data$sample.qc.missing.rate[match(related$ID2, data$eid)]
#   related$max_miss <- pmax(related$miss1, related$miss2)
# 
#   # Remove according to rule above
#   related$id_remove <- ifelse(is.na(related$miss1) & is.na(related$miss2),
#                               related$ID2,
#                               ifelse(is.na(related$miss1), related$ID1,
#                                      ifelse(is.na(related$miss2), related$ID2,
#                                             ifelse(related$miss1 ==
#                                                      related$max_miss, related$ID1,
#                                                    ifelse(related$miss2 ==
#                                                             related$max_miss,
#                                                           related$ID2, "error")))))
# 
#   cleaned <- subset(data, !(data$eid %in% related$id_remove))
# 
#   sink(qc_log_file, append = T)
#   cat(paste("**FILTER** Relatedness pairs with errors: ",
#             length(which(related$id_remove == "error")), "\n",
#             "**FILTER** Individuals excluded because of relatedness: ",
#             nrow(data[data$eid %in% related$id_remove, ]), "\n",
#             "REMAINING NOT RELATED: ", nrow(cleaned), "\n\n", sep = ""))
#   sink()
# 
#   return(cleaned)
# 
# }
# 
# qc_kinship_table <- function(data, qc_log_file) {
#   
#   cleaned <- subset(data, !is.na(data$excluded.from.kinship.inference) & 
#                       data$excluded.from.kinship.inference == 0)
#   
#   sink(qc_log_file, append = T)
#   cat(paste("**FILTER** Excluded from kinship inference: ", 
#             nrow(data) - nrow(cleaned), 
#             " ; REMAINING: ", nrow(cleaned), "\n", sep = ""))
#   sink()
#   
#   return(cleaned)
#   
# }


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

qcd_slopes <- lapply(PHENOTYPES, function (p) {
  res <- lapply(GWAS_STRATA, function (s) {
    # Stratum data
    data <- for_gen_QC[[p]][[s]]
    qc_log_file <- paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/log_files/genotyping_QC/",
                          p, "_", s, ".txt")
    
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
    # cleaned <- qc_excess_related(cleaned, qc_log_file)
    # cleaned <- qc_related(cleaned, qc_log_file)
    # cleaned <- qc_kinship_table(cleaned, qc_log_file)
    cleaned <- ukb_recommended_excl(cleaned, qc_log_file)
    
    return (cleaned)
  })
  names(res) <- GWAS_STRATA
  return (res)
})
names(qcd_slopes) <- PHENOTYPES

# Print list of all individuals that pass QC 

passed_QC <- lapply(PHENOTYPES, function (p) {
  eids <- lapply(GWAS_STRATA, function (s) {
    return (qcd_slopes[[p]][[s]]$eid)
  })
  return (unique(unlist(eids)))
})
passed_QC <- unique(unlist(passed_QC))

passed_QC <- data.frame("FID" = passed_QC, "IID" = passed_QC)
write.table(passed_QC, "GWAS/eids_passed_QC.txt", 
            sep = "\t", col.names = F, row.names = F, quote = F)

# Split into full population and gainers ----

split_slopes <- lapply(PHENOTYPES, function (p) {
  res <- lapply(GWAS_STRATA, function (s) {
    df <- qcd_slopes[[p]][[s]]
    df_list <- list()
    df_list$all <- df
    df_list$gainers <- df[df$gainer, ]
    return (df_list)
  })
  names(res) <- GWAS_STRATA
  return (res)
})
names(split_slopes) <- PHENOTYPES
GROUPS <- c("all", "gainers")

# RINT residuals and gather final phenotypes ----

RINTed <- lapply(PHENOTYPES, function (p) {
  res <- lapply(GWAS_STRATA, function (s) {
    res <- lapply(GROUPS, function (gp) {
      df <- split_slopes[[p]][[s]][[gp]][, c("eid", "FUyrs", "residual", 
                                             "genotyping.array")]
      # Calculate FUyr quartiles
      df$FUyr_quartile <- cut(df$FUyrs, quantile(df$FUyrs), 
                              include.lowest = T, labels = 1:4)
      # Calculate Z-score within each FUyr quartile
      df <- df %>% group_by(FUyr_quartile) %>% 
        mutate(Zscore = (residual - mean(residual)) / sd(residual))
      # Apply RINT on residuals
      df$RINT_resids <- qnorm((rank(df$residual) - 0.5) / 
                                sum(!is.na(df$residual)))
      # Apply RINT on Z-scores (normalised residuals)
      df$RINT_zscores <- qnorm((rank(df$Zscore) - 0.5) / sum(!is.na(df$Zscore)))
      
      # SELECTED COLUMNS FOR GWAS - IF STRATEGY CHANGES, REDO THESE
      df <- df[, c("eid", "eid", "genotyping.array", "FUyr_quartile",
                   "RINT_resids")]
      
      colnames(df) <- c("FID", "IID", "genotyping_array", "FUyr_quartile",
                        "RINTed_resids")
      
      write.table(df, paste0("GWAS/slope_files/", p, "_", s, "_", gp, ".txt"),
                  sep = "\t", quote = F, row.names = F)
      
      return (df)
    })
    names(res) <- GROUPS
    return (res)
  })
  names(res) <- GWAS_STRATA
  return (res)
})
names(RINTed) <- PHENOTYPES

saveRDS(RINTed, "/well/lindgren/UKBIOBANK/samvida/adiposity/RINTed_GWAS_phenotypes.rds")