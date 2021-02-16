# Author: Samvida S. Venkatesh
# Date: 10/02/21

library(tidyverse)
library(lubridate)

# Read data ----

# Read annotated gp_clinical file
gp_clinical <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/gp_clinical_annotated.txt",
                          sep = "\t", header = T, comment.char = "$",
                          stringsAsFactors = F)

# Remove individuals who have withdrawn consent
withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20210201.csv", 
                        header = F)
gp_clinical <- subset(gp_clinical, !gp_clinical$eid %in% withdrawn$V1)

# Read main UKB phenotypes file
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

# Extract adiposity read codes from primary care ----

PHENOTYPES <- c("weight", "waist", "BMI", "WHR")
V2_CODES <- c("22A..", "22N0.", "22K..", "22N7.")
V3_CODES <- c("22A..", "Xa041", "22K..", "X76CX")

names(V2_CODES) <- PHENOTYPES
names(V3_CODES) <- PHENOTYPES

adiposity_gp <- lapply(PHENOTYPES, function (p) {
  # Subset codes of interest
  df <- subset(gp_clinical, gp_clinical$read_2 == V2_CODES[p] |
                 gp_clinical$read_3 == V3_CODES[p])
  # Retain only the columns of interest
  # Only retain value1 as value2 and value3 sometimes have misleading values,
  # ex. BMI instead of weight rather than taking the mean or median
  df <- df[, c("eid", "data_provider", "event_dt", "age_years", "value1")]
  colnames(df) <- c("eid", "data_provider", "event_dt", "age_event", "value")
  
  df$data_provider <- paste0("GP", df$data_provider)
  # Convert all columns to standard type
  df$eid <- factor(df$eid)
  df$data_provider <- factor(df$data_provider)
  df$event_dt <- as.Date(df$event_dt, "%Y-%m-%d")
  df$age_event <- as.numeric(df$age_event)
  df$value <- as.numeric(df$value)
  return (df)
})
names(adiposity_gp) <- PHENOTYPES

# Extract adiposity field codes from UKB ----

FIELD_CODES <- list("weight" = c("f.21002.0.0", "f.21002.1.0", "f.21002.2.0"),
                    "waist" = c("f.48.0.0", "f.48.1.0", "f.48.2.0"),
                    "BMI" = c("f.21001.0.0", "f.21001.1.0", "f.21001.2.0"),
                    "WHR" = c("f.48.0.0", "f.48.1.0", "f.48.2.0",
                              "f.49.0.0", "f.49.1.0", "f.49.2.0"))
EVENT_DATE_CODES <- c("f.53.0.0", "f.53.1.0", "f.53.2.0")
AGE_CODES <- c("f.21003.0.0", "f.21003.1.0", "f.21003.2.0")
ASSMT_CENTRE_CODES <- c("f.54.0.0", "f.54.1.0", "f.54.2.0")

adiposity_ukb <- lapply(PHENOTYPES, function (p) {
  # Subset codes of interest
  # Also include event date, age at measurement, and assessment centre
  if (p == "WHR") {
    df <- pheno[, c("f.eid", EVENT_DATE_CODES, AGE_CODES,
                    ASSMT_CENTRE_CODES,
                    FIELD_CODES[[p]])]
    colnames(df) <- c("eid", 
                      "event_dt.1", "event_dt.2", "event_dt.3",
                      "age_event.1", "age_event.2", "age_event.3",
                      "data_provider.1", "data_provider.2", "data_provider.3",  
                      "waist.1", "waist.2", "waist.3", 
                      "hip.1", "hip.2", "hip.3")
    df <- pivot_longer(df, cols = -1, names_to = c(".value", "set"),
                       names_pattern = "(.+).(.+)", values_drop_na = T) %>%
      mutate(value = waist/hip)
    df$data_provider <- paste0("UKBB", df$data_provider) 
    df <- df[, c("eid", "data_provider", "event_dt", "age_event", "value")]
    # Convert all columns to standard type
    df$eid <- factor(df$eid)
    df$data_provider <- factor(df$data_provider)
    df$event_dt <- as.Date(df$event_dt, "%Y-%m-%d")
    df$age_event <- as.numeric(df$age_event)
    df$value <- as.numeric(df$value)
  } else {
    df <- pheno[, c("f.eid", EVENT_DATE_CODES, AGE_CODES,
                    ASSMT_CENTRE_CODES,
                    FIELD_CODES[[p]])]
    # Match column names with GP data
    colnames(df) <- c("eid", 
                      "event_dt.1", "event_dt.2", "event_dt.3",
                      "age_event.1", "age_event.2", "age_event.3",
                      "data_provider.1", "data_provider.2", "data_provider.3",  
                      "value.1", "value.2", "value.3")
    df <- pivot_longer(df, cols = -1, names_to = c(".value", "set"),
                       names_pattern = "(.+).(.+)", values_drop_na = T)
    df$data_provider <- paste0("UKBB", df$data_provider) 
    df <- df[, c("eid", "data_provider", "event_dt", "age_event", "value")]
    # Convert all columns to standard type
    df$eid <- factor(df$eid)
    df$data_provider <- factor(df$data_provider)
    df$event_dt <- as.Date(df$event_dt, "%Y-%m-%d")
    df$age_event <- as.numeric(df$age_event)
    df$value <- as.numeric(df$value)
  }
  return (df)
})
names(adiposity_ukb) <- PHENOTYPES

# Combine GP and UKB information ----

adiposity <- lapply(PHENOTYPES, function (p) {
  # Combine GP and UKB data
  df <- bind_rows(adiposity_gp[[p]], adiposity_ukb[[p]]) %>% 
    arrange(eid, event_dt)
  # Drop any NA values and only keep individuals with > 1 measurement
  df <- df[complete.cases(df$age_event, df$value), ]
  keep_ids <- df %>% group_by(eid) %>% count()
  keep_ids <- keep_ids$eid[keep_ids$n > 1]
  df <- subset(df, df$eid %in% keep_ids)
  return (df)
})
names(adiposity) <- PHENOTYPES

saveRDS(adiposity, "/well/lindgren/UKBIOBANK/samvida/adiposity/raw_adiposity.rds")

# Extract demographic and covariate information ----

# Get IDs of individuals on whom we have repeat adiposity data
ADIPO_EIDS <- unique(unlist(lapply(adiposity, function (df) df$eid )))

# UKB covariates: sex, ancestry, and median height
# Sex: 0 = F, 1 = M
# Ancestry: 1xxx = White, 2xxx = Mixed, 3XXX = Asian, 4XXX = Black,
# 5XXX = Chinese, 6XXX = Other
# Height  

covars <- pheno[pheno$f.eid %in% ADIPO_EIDS, 
                c("f.eid", "f.31.0.0", 
                  "f.21000.0.0", "f.21000.1.0", "f.21000.2.0",
                  "f.50.0.0", "f.50.1.0", "f.50.2.0")]
colnames(covars) <- c("eid", "sex",
                      "ancestry.1", "ancestry.2", "ancestry.3",
                      "height.1", "height.2", "height.3")

covars$sex <- ifelse(covars$sex == 0, "F", "M")

# Write function to name ancestry from codes
nameAncestry <- function (x) {
  ifelse(all(grepl("^1", x[!is.na(x)])), "white", 
         ifelse(all(grepl("^1", x[!is.na(x)])), "mixed",
                ifelse(all(grepl("^1", x[!is.na(x)])), "asian",
                       ifelse(all(grepl("^1", x[!is.na(x)])), "black",
                              ifelse(all(grepl("^1", x[!is.na(x)])), "chinese",
                                     ifelse(all(grepl("^1", x[!is.na(x)])), 
                                            "other",
                                            "missing or inconsistent"))))))
}
covars$ancestry <- apply(select(covars, starts_with("ancestry")), 1, 
                         function (x) nameAncestry(x) )

# Calculate median height
covars$height <- apply(select(covars, starts_with("height")), 1, 
                       function (x) median(as.numeric(x), na.rm = T) )

covars <- covars[, c("eid", "sex", "ancestry", "height")]

write.table(covars, "/well/lindgren/UKBIOBANK/samvida/adiposity/raw_adiposity_covars.txt",
            sep = "\t", quote = F, row.names = F)
