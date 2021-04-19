# Author: Samvida S. Venkatesh
# Date: 20/02/21
# ADAPTED FROM ADIPOSITY CHANGE CONSORTIUM SOP

library(lme4)
library(tidyverse)
#theme_set(theme_bw())

#set.seed(200221)

# Read files ----

adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/visually_QCd_adiposity.rds")
covars <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/model_covariates.rds")
PHENOTYPES <- names(adiposity)

NPCs <- 21
PCs <- paste0("PC", 1:NPCs)

STRATA <- lapply(covars, function (x) { x %>% distinct(ancestry, sex) }) 

# Create individual-level raw slopes ----

slope_summaries <- lapply(PHENOTYPES, function (p) {
  sub_covars <- covars[[p]]
  adipo <- adiposity[[p]]
  
  # One last check to remove any individuals without multiple measures
  keep_eids <- sub_covars$eid[sub_covars$FU_n > 1]
  adipo <- subset(adipo, adipo$eid %in% keep_eids)
  
  slopes <- lapply(1:nrow(STRATA[[p]]), function (si) {
    # Subset adiposity data relevant to stratum
    si_eids <- sub_covars$eid[sub_covars$sex == STRATA[[p]]$sex[si] &
                                sub_covars$ancestry == STRATA[[p]]$ancestry[si]]
    dat <- subset(adipo, adipo$eid %in% si_eids)
    # mixed effects model for adiposity on age
    mixed_model <- lmer(value ~ age_event + (age_event | eid), 
                        data = dat, REML = F)
    # produce a new variable for the random slope
    rs <- ranef(mixed_model)$eid$age_event
    # create a new variable for the fixed effect of age on adiposity
    bl <- fixef(mixed_model)[2]
    # sum the random and fixed effect to get individual slope (BLUP)
    slope <- rs + bl
    # create slopes to add to covariate summary data frame
    slopes <- data.frame(eid = rownames(ranef(mixed_model)$eid),
                         raw_slope = slope)
    
    # Flag outlier slopes > 5 S.D. away from mean
    # Remove values +/- 5 S.D. away from the strata mean slope
    popn_mean <- mean(slopes$raw_slope, na.rm = T)
    popn_sd <- sd(slopes$raw_slope, na.rm = T)
    max_outlier <- popn_mean + 5*popn_sd
    min_outlier <- popn_mean - 5*popn_sd
    
    slopes$slope_outlier_flag <- slopes$raw_slope > max_outlier |
      slopes$raw_slope < min_outlier
    
    # Report QC metrics
    sink(paste0("log_files/raw_slopes_QC_", p, ".txt"), append = T)
    cat(paste0("Strata: ", si, " ", 
               STRATA[[p]]$ancestry[si], " ", STRATA[[p]]$sex[si], "\n",
               "**FILTER** EXCLUDED, Raw slope > 5 S.D. away from strata mean: ", 
               sum(slopes$slope_outlier_flag), "\n"))
    sink()
    
    slopes$eid <- as.character(slopes$eid)
    
    return (slopes)
    
  })
  
  # Bind all strata slopes together
  slopes <- bind_rows(slopes)
  # Add to individual-level covariates data
  res <- merge(sub_covars, slopes, by = "eid")
  
  # Exclude outlier-flagged slopes
  res <- subset(res, !res$slope_outlier_flag)
  
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
names(slope_summaries) <- PHENOTYPES

# Save raw slopes list with sex-combined dataframes
saveRDS(slope_summaries, "/well/lindgren/UKBIOBANK/samvida/adiposity/raw_slopes_and_covars.rds")

# Stratify adiposity data and save QCd data ----

QCd_adipo <- lapply(PHENOTYPES, function (p) {
  adipo_sub <- adiposity[[p]]
  res <- lapply(slope_summaries[[p]], function (strata) {
    return (adipo_sub[adipo_sub$eid %in% strata$eid, ])
  })
  names(res) <- names(slope_summaries[[p]])
  return (res)
})
names(QCd_adipo) <- PHENOTYPES

saveRDS(QCd_adipo, "/well/lindgren/UKBIOBANK/samvida/adiposity/stratified_adiposity.rds")

# # Plot trajectories of individuals with extreme slopes ----
# 
# PHENOTYPE_LABELS <- c("weight (kg)", 
#                       "waist circumference (cm)",
#                       "BMI (kg/m^2)",
#                       "waist-hip ratio")
# names(PHENOTYPE_LABELS) <- PHENOTYPES
# 
# plot_extremes <- lapply(PHENOTYPES, function (p) {
#   strata_names <- names(slope_summaries[[p]])
#   
#   res <- lapply(strata_names, function (s) {
#     # Plot trajectories, randomly sampling from population in strata
#     # and also randomly sampling from IDs with extreme slopes
#     covars_df <- slope_summaries[[p]][[s]]
#     
#     RAND_non_outlier_IDs <- 
#       sample(covars_df$eid[!covars_df$slope_outlier_flag], 
#              min(20, length(covars_df$eid[!covars_df$slope_outlier_flag])), 
#              replace = F)
#     RAND_outlier_IDs <- 
#       sample(covars_df$eid[covars_df$slope_outlier_flag], 
#              min(20, length(covars_df$eid[covars_df$slope_outlier_flag])), 
#              replace = F)
#     plot_df <- subset(adiposity[[p]], 
#                       adiposity[[p]]$eid %in% c(RAND_non_outlier_IDs,
#                                                 RAND_outlier_IDs)) %>%
#       mutate(outlier = eid %in% RAND_outlier_IDs)
#     
#     res <- ggplot(plot_df, aes(x = age_event, y = value, group = eid,
#                                color = outlier)) +
#       geom_line() +
#       geom_point() +
#       scale_color_manual(values = c("#ebebeb", "#000000"), guide = F) +
#       labs(x = "age (years)", y = PHENOTYPE_LABELS[p], title = s)
#     return (res)
#   })
#   names(res) <- strata_names
#   pdf(paste0("plots/QC/extreme_slopes_", p, ".pdf"), onefile = T)
#   print(res)
#   dev.off()
#   return (res)
# })
# names(plot_extremes) <- PHENOTYPES
# 
# # Plot all adiposity trajectories for individuals with extreme slopes ----
# 
# # Get EIDs of individuals who are flagged as extreme outliers in each 
# # adiposity trait 
# 
# extreme_eids <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(slope_summaries[[p]], function (slopes) {
#     unique(slopes$eid[slopes$slope_outlier_flag])
#   })
#   res <- unlist(res)
#   return (res)
# })
# names(extreme_eids) <- PHENOTYPES
# 
# pdf("plots/QC/multivariate_extreme_slopes.pdf", onefile = T)
# upset(fromList(extreme_eids), order.by = "freq")
# dev.off()
# 
# # Get EIDs of all individuals
# all_eids <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(slope_summaries[[p]], function (slopes) {
#     unique(slopes$eid)
#   })
#   res <- unlist(res)
#   return (res)
# })
# names(all_eids) <- PHENOTYPES
# # Get EIDs for individuals for whom we have all phenotypes measured
# all_traits_eids <- Reduce(intersect, all_eids)
# 
# # Subset individuals who are outliers in only weight
# weight_eids <- extreme_eids$weight[!extreme_eids$weight %in% 
#                                      unlist(extreme_eids[2:4])]
# weight_eids <- intersect(all_traits_eids, weight_eids)
# # Subset individuals who are outliers only in BMI 
# BMI_eids <- extreme_eids$BMI[!extreme_eids$BMI %in% 
#                                      unlist(extreme_eids[c(1, 2, 4)])]
# BMI_eids <- intersect(all_traits_eids, BMI_eids)
# # Subset individuals who are outliers in both weight and BMI
# both_eids <- intersect(extreme_eids$weight, extreme_eids$BMI)
# both_eids <- intersect(all_traits_eids, both_eids)
# 
# # Plot trajectories for a sample of these individuals
# plot_eids <- c(sample(weight_eids, 1, replace = F),
#                sample(BMI_eids, 1, replace = F),
#                sample(both_eids, 1, replace = F))
# plot_df <- lapply(PHENOTYPES, function (p) {
#   res <- adiposity[[p]][adiposity[[p]]$eid %in% plot_eids, ]
#   res$trait <- p
#   return (res)
# })
# plot_df <- bind_rows(plot_df)
# 
# pdf("plots/QC/multivariate_trajectories_outlier_sample.pdf", onefile = T)
# ggplot(plot_df, aes(x = age_event, y = value, 
#                     group = eid, color = eid)) +
#   facet_wrap(~trait, nrow = 2, scales = "free_y") +
#   geom_line() +
#   geom_point() +
#   scale_color_brewer(palette = "Dark2", guide = F) +
#   labs(x = "age (years)", y = "Adiposity trait value")
# dev.off()
# 
# 
