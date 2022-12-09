# Author: Samvida S. Venkatesh
# Date: 19/04/2021

library(tidyverse)

# Read files ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/visually_QCd_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/model_covariates.rds")
PHENOTYPES <- names(adiposity)

NPCs <- 21
PCs <- paste0("PC", 1:NPCs)

STRATA <- lapply(covars, function (x) { x %>% distinct(ancestry, sex) }) 

# Create median trait measure ----

median_traits <- lapply(PHENOTYPES, function (p) {
  sub_covars <- covars[[p]]
  adipo <- adiposity[[p]]
  
  # One last check to remove any individuals without multiple measures
  keep_eids <- sub_covars$eid[sub_covars$FU_n > 1]
  adipo <- subset(adipo, adipo$eid %in% keep_eids)
  
  medians <- lapply(1:nrow(STRATA[[p]]), function (si) {
    # Subset adiposity data relevant to stratum
    si_eids <- sub_covars$eid[sub_covars$sex == STRATA[[p]]$sex[si] &
                                sub_covars$ancestry == STRATA[[p]]$ancestry[si]]
    dat <- subset(adipo, adipo$eid %in% si_eids)
    
    # Adiposity trait and age at midpoint of all measurements
    meds <- dat %>% group_by(eid) %>% 
      summarise(median_trait = value[ceiling(n()/2)],
                median_age = age_event[ceiling(n()/2)])
    
    # Flag outlier medians > 5 S.D. away from mean
    # Remove values +/- 5 S.D. away from the strata mean median
    popn_mean <- mean(meds$median_trait, na.rm = T)
    popn_sd <- sd(meds$median_trait, na.rm = T)
    max_outlier <- popn_mean + 5*popn_sd
    min_outlier <- popn_mean - 5*popn_sd
    
    meds$median_outlier_flag <- meds$median_trait > max_outlier |
      meds$median_trait < min_outlier
    
    # Report QC metrics
    sink(paste0("log_files/median_QC", p, ".txt"), append = T)
    cat(paste0("Strata: ", si, " ", 
               STRATA[[p]]$ancestry[si], " ", STRATA[[p]]$sex[si], "\n",
               "**FILTER** EXCLUDED, Median trait > 5 S.D. away from strata mean: ", 
               sum(meds$median_outlier_flag), "\n"))
    sink()
    
    meds$eid <- as.character(meds$eid)
    
    return (meds)
    
  })
  
  # Bind all strata medians together
  medians <- bind_rows(medians)
  # Add to individual-level covariates data
  res <- merge(sub_covars, medians, by = "eid")
  
  # Exclude outlier-flagged slopes
  res <- subset(res, !res$median_outlier_flag)
  
  # Split into sex- and ancestry-defined strata
  ss <- res %>% group_by(sex, ancestry) %>% group_split(.keep = T)
  names(ss) <- lapply(ss, function (x) paste(unique(x$ancestry), unique(x$sex),
                                             sep = "_"))
  # Also create a sex-combined stratum
  sc <- res %>% group_by(ancestry) %>% group_split(.keep = T)
  names(sc) <- lapply(sc, function (x) paste(unique(x$ancestry), "sexcomb",
                                             sep = "_"))
  res <- c(ss, sc)
  return (res)
})
names(median_traits) <- PHENOTYPES

# Save median traits list with sex-combined dataframes
saveRDS(median_traits, "/well/lindgren/UKBIOBANK/samvida/adiposity/covars_with_median_trait.rds")