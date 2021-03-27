# Author: Samvida S. Venkatesh
# Date: 16/02/21

library(tidyverse)
theme_set(theme_bw())
library(UpSetR)
library(RColorBrewer)

set.seed(160221)

# Read data ----
adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/QCd_adiposity.rds")

# covars <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/raw_adiposity_covars.txt",
#                      sep = "\t", header = T, comment.char = "$",
#                      stringsAsFactors = F)

PHENOTYPES <- names(adiposity)
PHENOTYPE_LABELS <- c("weight (kg)", 
                      "waist circumference (cm)",
                      "BMI (kg/m^2)",
                      "waist-hip ratio")
names(PHENOTYPE_LABELS) <- PHENOTYPES

# covars <- subset(covars, covars$eid %in% unique(unlist(lapply(adiposity, 
#                                                        function (df) df$eid ))))

# # Sanity check - age distribution of pregnancy flags ----
# 
# # There are mistakes in the GP data that lead to pregnancy flags in women
# # over the age of 60; remove these
# adiposity <- lapply(adiposity, function (df) {
#   df$preg_flag[df$age_event > 60] <- F
#   return (df)
# })
# 
# plot_preg_age_distribution <- lapply(PHENOTYPES, function (p) {
#   df <- adiposity[[p]]
#   preg_ages <- data.frame(preg_age = df$age_event[df$preg_flag])
#   res <- ggplot(preg_ages, aes(x = preg_age)) +
#     geom_density() +
#     labs(x = "age at pregnancy flag (years)", title = paste(p, "measured"), 
#          y = "density")
#   return (res)
# })
# 
# pdf("plots/QC/age_distribution_preg_flags.pdf", onefile = T)
# print(plot_preg_age_distribution)
# dev.off()
# 
# # Sanity check - adiposity in individuals with bariatric surgery flags ----
# 
# plot_bar_adipo_distribution <- lapply(PHENOTYPES, function (p) {
#   df <- adiposity[[p]]
#   EVER_BAR <- unique(df$eid[df$post_bariatric_surgery_flag])
#   bar_adipo <- data.frame(adipo = df$value[df$eid %in% EVER_BAR])
#   not_bar_adipo <- data.frame(adipo = df$value[!df$eid %in% EVER_BAR])
#   res <- ggplot(not_bar_adipo, aes(x = adipo)) +
#     geom_density(color = "#ebebeb", fill = "#ebebeb", 
#                  alpha = 0.4) +
#     geom_density(data = bar_adipo, aes(x = adipo),
#                  color = "#D95F02", fill = "#D95F02",
#                  alpha = 0.2) +
#     xlim(c(min(min(bar_adipo$adipo), min(not_bar_adipo$adipo)),
#            max(max(bar_adipo$adipo), max(not_bar_adipo$adipo)))) +
#     labs(x = PHENOTYPE_LABELS[p], y = "density")
#   return (res)
# })
# 
# pdf("plots/QC/adiposity_distribution_bariatric_surgery.pdf", onefile = T)
# print(plot_bar_adipo_distribution)
# dev.off()
# 
# # Plots ----
# 
# # Plot trajectories, randomly sampling from population (sex-separate)
# # And also randomly sampling from bariatric surgery and pregnancy-flagged EIDs
# 
# RANDOM_F <- sample(covars$eid[covars$sex == "F"], 15, replace = F)
# RANDOM_M <- sample(covars$eid[covars$sex == "M"], 15, replace = F)
# 
# plot_dfs <- lapply(adiposity, function (df) {
#   # Get flagged individuals
#   bar_flag <- as.character(unique(df$eid[which(df$post_bariatric_surgery_flag)]))
#   bar_flag <- data.frame(eid = bar_flag, 
#                          sex = covars$sex[match(bar_flag, covars$eid)])
#   preg_flag <- as.character(unique(df$eid[which(df$preg_flag)]))
#   # Randomly sample for plots
#   bar_random_F <- sample(bar_flag$eid[bar_flag$sex == "F"], 
#                          min(10, length(bar_flag$eid[bar_flag$sex == "F"])), 
#                          replace = F)
#   bar_random_M <- sample(bar_flag$eid[bar_flag$sex == "M"], 
#                          min(10, length(bar_flag$eid[bar_flag$sex == "M"])), 
#                          replace = F)
#   preg_random <- sample(preg_flag, 20, replace = F)
#   # Subset relevant IDs from dataframe
#   plot_df_F <- subset(df, df$eid %in% 
#                         c(RANDOM_F, bar_random_F, preg_random)) %>%
#     mutate(sex = "F") 
#   plot_df_M <- subset(df, df$eid %in% 
#                         c(RANDOM_M, bar_random_M)) %>% mutate(sex = "M")
#   plot_df <- bind_rows(plot_df_F, plot_df_M)
#   
#   plot_df$EVER_BAR <- plot_df$eid %in% bar_flag$eid
#   plot_df$EVER_PREG <- plot_df$eid %in% preg_flag
#   return (plot_df)
# })
# 
# ## Bariatric surgery plots ----
# 
# bariatric_plots <- lapply(PHENOTYPES, function (p) {
#   df <- plot_dfs[[p]]
#   res <- ggplot(subset(df, !df$EVER_BAR), 
#                 aes(x = age_event, y = value, group = eid)) +
#     facet_wrap(~sex, nrow = 2) +
#     # plot non-surgery samples
#     geom_line(color = "#ebebeb") +
#     geom_point(color = "#ebebeb") +
#     # plot surgery samples
#     geom_line(data = subset(df, df$EVER_BAR), 
#               aes(color = post_bariatric_surgery_flag)) +
#     geom_point(data = subset(df, df$EVER_BAR), 
#                aes(color = post_bariatric_surgery_flag)) +
#     scale_color_manual(values = c("#000000", "#D95F02"), guide = F) +
#     ylim(min(df$value), max(df$value)) +
#     labs(x = "age (years)", y = PHENOTYPE_LABELS[p])
#   return (res)
# })
# 
# pdf("plots/QC/bariatric_surgery_trajectories.pdf", onefile = T)
# print(bariatric_plots)
# dev.off()
# 
# ## Pregnancy plots ----
# 
# preg_plots <- lapply(PHENOTYPES, function (p) {
#   df <- plot_dfs[[p]]
#   df <- subset(df, df$sex == "F" & !df$EVER_BAR)
#   
#   res <- ggplot(subset(df, !df$EVER_PREG), 
#                 aes(x = age_event, y = value, group = eid)) +
#     # plot never pregnant samples
#     geom_line(color = "#ebebeb") +
#     geom_point(color = "#ebebeb") +
#     # plot pregnant samples
#     geom_line(data = subset(df, df$EVER_PREG), aes(color = preg_flag)) +
#     geom_point(data = subset(df, df$EVER_PREG), aes(color = preg_flag)) +
#     scale_color_manual(values = c("#000000", "#D95F02"), guide = F) +
#     labs(x = "age (years)", y = PHENOTYPE_LABELS[p])
#   return (res)
# })
# 
# pdf("plots/QC/pregnancy_trajectories.pdf", onefile = T)
# print(preg_plots)
# dev.off()

# Find and remove large jumps ----

## Calculate jumps ----

cleaned <- adiposity

# Calculate jump as time-adjusted fold-change
jumps <- lapply(cleaned, function (df) {
  
  # If two measurements are taken at the same age, 
  # remove the one farther from the overall median value
  dups <- df %>% group_by(eid) %>%
    mutate(dist_to_med = abs(value - median(value)),
           # check if there are two values at the same age
           dup_marker = duplicated(age_event) | duplicated(age_event, fromLast = T))
  
  res <- subset(dups, !dups$dup_marker)
  if (any(dups$dup_marker)) {
    no_dups <- subset(dups, dups$dup_marker) %>% group_by(eid, age_event) %>%
      mutate(remove = dist_to_med == max(dist_to_med))
    no_dups <- subset(no_dups, !no_dups$remove)
    res <- bind_rows(res, no_dups)
  }
  
  res <- res %>% group_by(eid) %>% 
    # Ensure ordered
    arrange(age_event, .by_group = T) %>% 
    mutate(FC = abs(value - lag(value)) / lag(value),
           age_change = age_event - lag(age_event),
           jump = ifelse(FC == 0, NA, log2(FC/age_change)))
  
  # Mark jump as extreme if > 5 S.D. larger than mean jump
  mean_jump <- mean(res$jump, na.rm = T) 
  sd_jump <- sd(res$jump, na.rm = T)
  res$extreme <- res$jump > (mean_jump + (5*sd_jump))
  
  return (res)
})

# ## Visually inspect extreme jumps ----

# Visually inspect random sample of trajectories with large jumps
jump_plots <- lapply(PHENOTYPES, function (p) {
  df <- jumps[[p]]
  random_nonextreme <- sample(unique(df$eid[!df$extreme]),
                              min(15, length(unique(df$eid[!df$extreme]))),
                              replace = F)
  random_extreme <- sample(unique(df$eid[df$extreme]),
                           min(15, length(unique(df$eid[df$extreme]))),
                           replace = F)

  res <- ggplot(subset(df, df$eid %in% random_nonextreme),
                aes(x = age_event, y = value, group = eid)) +
    # plot non-extreme samples
    geom_line(color = "#ebebeb") +
    geom_point(color = "#ebebeb") +
    # plot extreme samples
    geom_line(data = subset(df, df$eid %in% random_extreme),
              color = "#000000") +
    geom_point(data = subset(df, df$eid %in% random_extreme),
               aes(color = extreme)) +
    scale_color_manual(values = c("#000000", "#D95F02"), guide = F) +
    xlim(c(min(df$age_event), max(df$age_event))) +
    ylim(c(min(df$value), max(df$value))) +
    labs(x = "age (years)", y = PHENOTYPE_LABELS[p])
  return (res)
})

pdf("plots/QC/extreme_jumps_round1.pdf", onefile = T)
print(jump_plots)
dev.off()

# ## Remove values causing extreme jumps ----

cleaned <- lapply(jumps, function (df) {
  
  if (any(df$extreme, na.rm = T)) {
    extremes <- df[df$eid %in% df$eid[df$extreme], ]
    # Remove the point farthest from the median (that causes the jump)
    extremes <- extremes %>% group_by(eid) %>%
      mutate(median_value = median(value),
             dist_to_med = abs(value - median_value),
             remove = dist_to_med == max(dist_to_med))
    extremes <- subset(extremes, !extremes$remove)
    cleaned <- bind_rows(df[!df$eid %in% df$eid[df$extreme], ],
                         extremes[, colnames(df)])
  } else cleaned <- df
  
  return (cleaned)
})

# ITERATE THE ABOVE STEPS UNTIL VISUAL INSPECTION YIELDS NO
# OBVIOUS DATA ENTRY ERRORS IN JUMPS

# Data metrics before and after cleaning ----

lapply(PHENOTYPES, function (p) {
  
  orig_obs <- dim(adiposity[[p]])[1]
  orig_indivs <- length(unique(adiposity[[p]]$eid))
  
  cleaned_obs <- dim(cleaned[[p]])[1]
  cleaned_indivs <- length(unique(cleaned[[p]]$eid))
  
  sink(paste0("log_files/visual_qc_log_file_", p, ".txt"), 
       append = T)
  cat(paste0("**FILTER** EXCLUDED, Unreasonable jump by visual inspection: ", 
             "\n",
             "\t", "Number of measurements = ", 
             orig_obs - cleaned_obs, "\n",
             "\t", "Number of individuals = ", 
             orig_indivs - cleaned_indivs, "\n"))
  sink()
  
})

cleaned <- lapply(cleaned, function (df) {
  df <- df[, c("eid", "data_provider", "event_dt",
               "age_event", "value")]
  df$eid <- as.character(df$eid)
  return (df)
})

saveRDS(cleaned, "visually_QCd_adiposity.rds")

# Multi-variate data availability after cleaning ----

eid_lists <- lapply(cleaned, function (x) { unique(x$eid) })
pdf("plots/multivariate/adiposity_frequency.pdf", onefile = T)
upset(fromList(eid_lists), order.by = "freq")
dev.off()

# Trajectories for individuals who have all 4 traits measured - randomly sample
# and plot 5 individuals

# Get EIDs for individuals for whom we have all phenotypes measured

all_traits_eids <- Reduce(intersect, eid_lists)

plot_eids <- sample(all_traits_eids, 5, replace = F)
plot_df <- lapply(PHENOTYPES, function (p) {
  res <- adiposity[[p]][adiposity[[p]]$eid %in% plot_eids, ]
  res$trait <- p
  return (res)
})
plot_df <- bind_rows(plot_df)

pdf("plots/multivariate/trajectories_sample.pdf", onefile = T)
ggplot(plot_df, aes(x = age_event, y = value, 
                    group = eid, color = eid)) +
  facet_wrap(~trait, nrow = 2, scales = "free_y") +
  geom_line() +
  geom_point() +
  scale_color_brewer(palette = "Dark2", guide = F) +
  labs(x = "age (years)", y = "Adiposity trait value")
dev.off()
