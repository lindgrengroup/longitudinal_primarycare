# Author: Samvida S. Venkatesh
# Date: 02/03/21

library(tidyverse)
#library(treemap)
library(RColorBrewer)
theme_set(theme_bw())

set.seed(020321)

resdir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/results/"

# Read adiposity and covariate data ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]
covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Construct demographic characteristics tables -----

demo_table <- lapply(PHENOTYPES, function (p) {
  # Count number of data points by sex and ancestry
  df <- dat[[p]] %>% left_join(general_covars[, c("eid", "sex")], 
                               by = "eid")
  obs_df <- df %>% 
    group_by(sex) %>% 
    summarise(n_obs = n())
  
  write.table(obs_df, paste0(resdir, "descriptive_factors/obs_counts_by_sex_", p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  
  # Covars data
  covars_df <- covars[[p]] %>% left_join(general_covars, by = "eid")
  summary_df <- covars_df %>% group_by(sex) %>% 
    summarise(count = n(), 
              mean_FU_n = mean(FU_n), sd_FU_n = sd(FU_n),
              to_print_FU_n = paste0(signif(mean_FU_n, 3), " (", signif(sd_FU_n, 3), ")"),
              mean_FUyrs = mean(FUyrs), sd_FUyrs = sd(FUyrs),
              to_print_FUyrs = paste0(signif(mean_FUyrs, 3), " (", signif(sd_FUyrs, 3), ")"),
              mean_bl_age = mean(baseline_age), 
              sd_bl_age = sd(baseline_age),
              to_print_bl_age = paste0(signif(mean_bl_age, 3), " (", signif(sd_bl_age, 3), ")"),
              median_bl_trait = median(baseline_trait), 
              iqr_bl_trait = paste(signif(quantile(baseline_trait, 0.25), 3), 
                                   signif(quantile(baseline_trait, 0.75), 3),
                                   sep = ", "),
              to_print_bl_trait = paste0(signif(median_bl_trait, 3), " (", iqr_bl_trait, ")")) %>%
    ungroup() %>% 
    mutate(percent = count / sum(count),
           to_print_n = paste0(count, " (", signif(percent*100, 3), "%)")) %>%
    select(all_of(c("sex", 
                    "to_print_n", "to_print_FU_n", "to_print_FUyrs", 
                    "to_print_bl_age", "to_print_bl_trait")))
  
  write.table(summary_df, paste0(resdir, "descriptive_factors/descriptive_factors_", p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  return (summary_df)
})

## Non-WB ancestry ----

dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/non_wb_gp_main_data_passed_longit_filter.rds")[PHENOTYPES]

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Construct demographic characteristics tables -----

demo_table <- lapply(PHENOTYPES, function (p) {
  # Count number of data points by sex and ancestry
  df <- dat[[p]] %>% left_join(general_covars[, c("eid", "sex", "ancestry")], 
                               by = "eid")
  obs_df <- df %>% 
    group_by(sex, ancestry) %>% 
    summarise(n_obs = n())
  
  write.table(obs_df, paste0(resdir, "descriptive_factors/non_wb_obs_counts_by_sex_", p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  
  # Covars data
  covars_df <- dat[[p]] %>% group_by(eid) %>% arrange(age_event) %>% 
    summarise(FU_n = n(), 
              FUyrs = last(age_event) - first(age_event),
              baseline_age = first(age_event),
              baseline_trait = first(value)) 

  covars_df <- covars_df %>% left_join(general_covars, by = "eid")
  summary_df <- covars_df %>% group_by(sex, ancestry) %>% 
    summarise(count = n(), 
              mean_FU_n = mean(FU_n), sd_FU_n = sd(FU_n),
              to_print_FU_n = paste0(signif(mean_FU_n, 3), " (", signif(sd_FU_n, 3), ")"),
              mean_FUyrs = mean(FUyrs), sd_FUyrs = sd(FUyrs),
              to_print_FUyrs = paste0(signif(mean_FUyrs, 3), " (", signif(sd_FUyrs, 3), ")"),
              mean_bl_age = mean(baseline_age), 
              sd_bl_age = sd(baseline_age),
              to_print_bl_age = paste0(signif(mean_bl_age, 3), " (", signif(sd_bl_age, 3), ")"),
              median_bl_trait = median(baseline_trait), 
              iqr_bl_trait = paste(signif(quantile(baseline_trait, 0.25), 3), 
                                   signif(quantile(baseline_trait, 0.75), 3),
                                   sep = ", "),
              to_print_bl_trait = paste0(signif(median_bl_trait, 3), " (", iqr_bl_trait, ")")) %>%
    ungroup() %>% group_by(ancestry) %>%
    mutate(percent = count / sum(count),
           to_print_n = paste0(count, " (", signif(percent*100, 3), "%)")) %>%
    select(all_of(c("sex", "ancestry",
                    "to_print_n", "to_print_FU_n", "to_print_FUyrs", 
                    "to_print_bl_age", "to_print_bl_trait")))
  
  write.table(summary_df, paste0(resdir, "descriptive_factors/non_wb_descriptive_factors_", p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  return (summary_df)
})

## THE REST OF THIS SCRIPT IS ARCHIVED -- PLOTS THAT WE DON'T USE IN THE MANUSCRIPT

# # Plot: ancestry and sex population counts ----
# 
# # demo_table <- lapply(PHENOTYPES, function (p) {
# #   res <- read.table(paste0("results/descriptive_factors/descriptive_factors_", p, ".txt"),
# #                     sep = "\t", header = T, stringsAsFactors = F)
# # })
# # names(demo_table) <- PHENOTYPES
# 
# popn_ancestry_treemap <- lapply(PHENOTYPES, function (p) {
#   df <- demo_table[[p]] %>% group_by(ancestry) %>% summarise(size = sum(count))
#   tiff(filename = paste0("plots/descriptive_factors/popn_ancestry_", p, ".tiff"),
#        height = 10, width = 10, units = "cm",
#        res = 600)
#   treemap(df, 
#           index = "ancestry", vSize = "size", type = "index",
#           palette = "Dark2", border.col = c("black"), border.lwds = 1,
#           title = "", fontsize.labels = 0)
#   dev.off()
#   return ()
# })
# 
# popn_sex_pies <- lapply(PHENOTYPES, function (p) {
#   df <- demo_table[[p]]
#   for (a in unique(df$ancestry)) {
#     res <- ggplot(subset(df, df$ancestry == a), 
#                   aes(x = "", y = count, fill = sex)) +
#       geom_bar(width = 1, stat = "identity") +
#       coord_polar("y", start = 0) +
#       scale_color_brewer(palette = "Dark2") +
#       theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
#             axis.title.y = element_blank(), panel.border = element_blank(), 
#             panel.grid = element_blank(), axis.ticks = element_blank(),
#             legend.position = "none")
#     tiff(filename = paste0("plots/descriptive_factors/popn_sex_", p, "_", 
#                            a, ".tiff"),
#          height = 10, width = 10, units = "cm",
#          res = 600)
#     print(res)
#     dev.off()
#   }
#   return ()
# })
# 
# # Box and violin plots of distributions ----
# 
# plot_distributions <- lapply(PHENOTYPES, function (p) {
#   pheno <- lapply(SS_STRATA, function (s) {
#     raw_slopes[[p]][[s]][, c("eid", "sex", "ancestry", "FU_n", "FUyrs", 
#                              "baseline_age", "height", 
#                              "baseline_BMI", "baseline_trait")]
#   })
#   df <- bind_rows(pheno)
#   
#   # Number of follow-up measures
#   FU_n <- ggplot(df, aes(x = ancestry, y = FU_n, fill = ancestry)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_boxplot(position = position_dodge(1)) +
#     scale_x_discrete(limits = 
#                        c("white", "asian", "other", "black", "mixed", "chinese")) + 
#     scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
#     scale_fill_brewer(palette = "Dark2") + 
#     labs(x = "ancestry", y = "# repeat measures") +
#     theme(legend.position = "none")
#   
#   # Number of years of follow up
#   FUyrs <- ggplot(df, aes(x = ancestry, y = FUyrs, fill = ancestry)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_boxplot(position = position_dodge(1)) +
#     scale_x_discrete(limits = 
#                        c("white", "asian", "other", "black", "mixed", "chinese")) + 
#     scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
#     scale_fill_brewer(palette = "Dark2") + 
#     labs(x = "ancestry", y = "follow-up length (years)") +
#     theme(legend.position = "none")
#   
#   # Baseline age
#   bl_age <- ggplot(df, aes(x = ancestry, y = baseline_age)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
#     geom_boxplot(width = 0.1) + 
#     scale_x_discrete(limits = 
#                        c("white", "asian", "other", "black", "mixed", "chinese")) + 
#     scale_y_continuous(limits = c(20, 80), breaks = seq(20, 80, by = 5)) +
#     scale_fill_brewer(palette = "Dark2") + 
#     labs(x = "ancestry", y = "baseline age (years)") +
#     theme(legend.position = "none")
#   
#   # Adult height
#   med_height <- ggplot(df, aes(x = ancestry, y = height)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
#     geom_boxplot(width = 0.1) + 
#     scale_x_discrete(limits = 
#                        c("white", "asian", "other", "black", "mixed", "chinese")) + 
#     scale_y_continuous(limits = c(110, 210), breaks = seq(110, 210, by = 10)) +
#     scale_fill_brewer(palette = "Dark2") + 
#     labs(x = "ancestry", y = "median adult height (cm)") +
#     theme(legend.position = "none")
#   
#   # Baseline BMI
#   bl_BMI <- ggplot(df, aes(x = ancestry, y = baseline_BMI)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
#     geom_boxplot(width = 0.1) + 
#     scale_x_discrete(limits = 
#                        c("white", "asian", "other", "black", "mixed", "chinese")) + 
#     scale_y_continuous(limits = c(10, 65), breaks = seq(10, 65, by = 5)) +
#     scale_fill_brewer(palette = "Dark2") + 
#     labs(x = "ancestry", y = "baseline BMI (kg/m2)") +
#     theme(legend.position = "none")
#   
#   # Baseline adiposity trait
#   bl_trait <- ggplot(df, aes(x = ancestry, y = baseline_trait)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
#     geom_boxplot(width = 0.1) + 
#     scale_x_discrete(limits = 
#                        c("white", "asian", "other", "black", "mixed", "chinese")) + 
#     scale_fill_brewer(palette = "Dark2") + 
#     labs(x = "ancestry", y = "baseline adiposity trait") +
#     theme(legend.position = "none")
#   
#   pdf(paste0("plots/descriptive_factors/covariate_distributions_", p, ".pdf"),
#       onefile = T)
#   print(FU_n)
#   print(FUyrs)
#   print(bl_age)
#   print(med_height)
#   print(bl_BMI)
#   print(bl_trait)
#   dev.off()
#   return ()
# })
# 
# # Summary plot of baseline age, FU length, FU measures -----
# 
# # demo_table <- lapply(PHENOTYPES, function (p) {
# #   res <- read.table(paste0("results/descriptive_factors/descriptive_factors_", 
# #                            p, ".txt"),
# #                     sep = "\t", header = T, stringsAsFactors = F)
# # })
# # names(demo_table) <- PHENOTYPES
# 
# demo_table <- lapply(PHENOTYPES, function (p) {
#   res <- demo_table[[p]]
#   res$adipo_trait <- p
#   return (res)
# })
# demo_table <- bind_rows(demo_table)
# # min-max normalise line thickness from 1:5
# demo_table$line_thick <- 1 + 
#   4*(demo_table$mean_FU_n - min(demo_table$mean_FU_n)) / 
#   (max(demo_table$mean_FU_n) - min(demo_table$mean_FU_n))
# 
# pdf("plots/descriptive_factors/summary_FU.pdf")
# ggplot(demo_table) +
#   facet_wrap(~ancestry+sex, ncol = 2) +
#   geom_segment(aes(x = mean_bl_age, xend = mean_bl_age + mean_FUyrs,
#                    y = adipo_trait, yend = adipo_trait, 
#                    size = mean_FU_n,
#                    color = ancestry)) +
#   scale_y_discrete(limits = c("WHR", "waist", "weight", "BMI"),
#                    labels = c("WHR", "WC", "weight", "BMI")) +
#   scale_x_continuous(limits = c(40, 65), breaks = seq(40, 65, by = 5)) +
#   scale_color_brewer(palette = "Dark2") + 
#   labs(x = "Age (years)", y = "Adiposity trait") +
#   theme(legend.position = "none", 
#         strip.background = element_blank(),
#         strip.text.x = element_blank())
# dev.off()
# 
# # Plot binned trajectories ----
# 
# # Calculate mean and SE in each 5-year interval in each strata
# 
# age_bin_cuts <- seq(20, 80, by = 5)
# 
# binned_adipo <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(SS_STRATA, function (s) {
#     df <- adiposity[[p]][[s]]
#     df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
#     res <- df %>% group_by(age_bin) %>% summarise(count = n(),
#                                                   mean_value = mean(value),
#                                                   se_value = sd(value)/sqrt(count))
#     res$ancestry <- strsplit(s, "_")[[1]][1]
#     res$sex <- strsplit(s, "_")[[1]][2]
#     return (res)
#   })
#   res <- bind_rows(res)
#   return (res)
# })
# names(binned_adipo) <- PHENOTYPES
# 
# traj_plots <- lapply(PHENOTYPES, function (p) {
#   res <- ggplot(binned_adipo[[p]], aes(x = age_bin, y = mean_value,
#                                        group = ancestry,
#                                        color = ancestry, fill = ancestry)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_point() +
#     geom_path() +
#     geom_ribbon(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
#                 alpha = 0.2) +
#     scale_fill_brewer(palette = "Dark2") +
#     scale_color_brewer(palette = "Dark2") +
#     labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
#          title = p)
#   return (res)
# })
# names(traj_plots) <- PHENOTYPES
# 
# pdf("plots/descriptive_factors/binned_mean_trajectories.pdf", onefile = T)
# print(traj_plots)
# dev.off()
# 
# # Plot sample trajectories ----
# 
# sample_trajectories <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(SS_STRATA, function (s) {
#     df <- adiposity[[p]][[s]]
#     ids <- sample(df$eid, 10, replace = F)
#     
#     plot_df <- long_adiposity[[p]]
#     plot_df <- subset(plot_df, plot_df$eid %in% ids)
#     
#     p <- ggplot(plot_df, aes(x = age_event, y = value, group = eid)) +
#       geom_line() +
#       geom_point() +
#       labs(x = "age (years)", y = "Adiposity trait value", 
#            title = s)
#     return (p)
#   })
#   names(res) <- SS_STRATA
#   pdf(paste0("plots/descriptive_factors/sample_trajectories_", p, ".pdf"),
#       onefile = T)
#   print(res)
#   dev.off()
#   return (res)
# })
# names(sample_trajectories) <- PHENOTYPES
# 
# # Plot change in trait vs mean trait (sample) ----
# 
# change_mean_dfs <- lapply(PHENOTYPES, function (pheno) {
#   res <- long_adiposity[[pheno]] %>% group_by(eid) %>%
#     arrange(age_event, .by_group = T) %>% 
#     mutate(mean_value_t1t2 = (value + lag(value))/2,
#            value_change = value - lag(value),
#            age_change = age_event - lag(age_event),
#            value_change_t1t2 = value_change / age_change)
# })
# 
# sample_trajectories <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(SS_STRATA, function (s) {
#     df <- adiposity[[p]][[s]]
#     ids <- sample(df$eid, 10, replace = F)
#     
#     plot_df <- change_mean_dfs[[p]]
#     plot_df <- subset(plot_df, plot_df$eid %in% ids)
#     
#     p <- ggplot(plot_df, aes(x = mean_value_t1t2, 
#                              y = value_change_t1t2, group = eid)) +
#       geom_line() +
#       geom_point() +
#       labs(x = "Mean adiposity trait value", y = "Change in adiposity / time", 
#            title = s)
#     return (p)
#   })
#   names(res) <- SS_STRATA
#   pdf(paste0("plots/descriptive_factors/sample_change_v_mean_", p, ".pdf"),
#       onefile = T)
#   print(res)
#   dev.off()
#   return (res)
# })
# names(sample_trajectories) <- PHENOTYPES
# 
# # Plot change in trait vs mean trait (binned) ----
# 
# binned_adipo <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(SS_STRATA, function (s) {
#     # Calculate mean and SE in each trait bin in each strata
#     df <- change_mean_dfs[[p]]
#     trait_bin_cuts <- seq(min(df$mean_value_t1t2, na.rm = T), 
#                           max(df$mean_value_t1t2, na.rm = T), 
#                           length.out = 11)
#     
#     df$trait_bin <- cut(df$mean_value_t1t2, trait_bin_cuts, include.lowest = T)
#     res <- df %>% group_by(trait_bin) %>% 
#       summarise(count = n(),
#                 mean_change_value = mean(value_change_t1t2, na.rm = T),
#                 se_value = sd(value_change_t1t2)/sqrt(count))
#     res$ancestry <- strsplit(s, "_")[[1]][1]
#     res$sex <- strsplit(s, "_")[[1]][2]
#     return (res)
#   })
#   res <- bind_rows(res)
#   return (res)
# })
# names(binned_adipo) <- PHENOTYPES
# 
# traj_plots <- lapply(PHENOTYPES, function (p) {
#   
#   res <- ggplot(binned_adipo[[p]], aes(x = age_bin, y = mean_value,
#                                        group = ancestry,
#                                        color = ancestry, fill = ancestry)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_point() +
#     geom_path() +
#     geom_ribbon(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
#                 alpha = 0.2) +
#     scale_fill_brewer(palette = "Dark2") +
#     scale_color_brewer(palette = "Dark2") +
#     labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
#          title = p)
#   return (res)
# })
# names(traj_plots) <- PHENOTYPES
# 
# pdf("plots/descriptive_factors/binned_mean_trajectories.pdf", onefile = T)
# print(traj_plots)
# dev.off()
# 
# # Construct raw slope summary table ----
# 
# rs_table <- lapply(PHENOTYPES, function (p) {
#   pheno <- lapply(SS_STRATA, function (s) {
#     raw_slopes[[p]][[s]][, c("eid", "sex", "ancestry", "raw_slope")]
#   })
#   df <- bind_rows(pheno)
#   summary_df <- df %>% group_by(sex, ancestry) %>% 
#     summarise(median_slope = median(raw_slope), 
#               iqr_slope = paste(quantile(raw_slope, 0.25), 
#                                 quantile(raw_slope, 0.75),
#                                 sep = ", "))
#   write.table(summary_df, paste0("results/raw_slopes/raw_slope_summaries_", p, ".txt"),
#               sep = "\t", quote = F, row.names = F)
#   return (summary_df)
# })
# 
# # Plot violin distribution of raw slopes ----
# 
# rs_distributions <- lapply(PHENOTYPES, function (p) {
#   pheno <- lapply(SS_STRATA, function (s) {
#     raw_slopes[[p]][[s]][, c("eid", "sex", "ancestry", "raw_slope")]
#   })
#   df <- bind_rows(pheno)
#   res <- ggplot(df, aes(x = ancestry, y = raw_slope)) +
#     facet_wrap(~sex, nrow = 2) +
#     geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
#     geom_boxplot(width = 0.1) + 
#     scale_x_discrete(limits = 
#                        c("white", "asian", "other", "black", "mixed", "chinese")) + 
#     scale_fill_brewer(palette = "Dark2") + 
#     labs(x = "ancestry", y = "raw slope", title = p) +
#     theme(legend.position = "none")
#   return (res)
# })
# pdf(paste0("plots/raw_slopes/raw_slopes_distributions.pdf"),
#     onefile = T)
# print(rs_distributions)
# dev.off()
# 
# # Plot trajectories of sample individuals at the tails of raw slope distribution ----
# 
# rs_tail_trajectories <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(SS_STRATA, function (s) {
#     rs_df <- raw_slopes[[p]][[s]]
#     top_ids <- rs_df$eid[rs_df$raw_slope > quantile(rs_df$raw_slope, 0.95)]
#     top_ids <- sample(top_ids, min(length(top_ids), 5), replace = F)
#     bottom_ids <- rs_df$eid[rs_df$raw_slope < quantile(rs_df$raw_slope, 0.05)]
#     bottom_ids <- sample(bottom_ids, min(length(bottom_ids), 5), replace = F)
#     
#     plot_df <- long_adiposity[[p]]
#     plot_df <- subset(plot_df, 
#                       plot_df$eid %in% top_ids | plot_df$eid %in% bottom_ids)
#     plot_df$tail <- ifelse(plot_df$eid %in% top_ids, "top", "bottom")
#     
#     p <- ggplot(plot_df, aes(x = age_event, y = value, 
#                              group = eid, color = tail)) +
#       geom_line() +
#       geom_point() +
#       scale_color_manual(values = c("top" = "#984EA3", "bottom" = "#E41A1C"), 
#                          guide = F) +
#       labs(x = "age (years)", y = "Adiposity trait value", 
#            title = s)
#     return (p)
#   })
#   names(res) <- SS_STRATA
#   pdf(paste0("plots/raw_slopes/tail_trajectories_", p, ".pdf"),
#       onefile = T)
#   print(res)
#   dev.off()
#   return (res)
# })
# names(rs_tail_trajectories) <- PHENOTYPES
# 
# # Plot mean (binned) trajectories by raw slope quartile ----
# 
# rs_quartile_trajectories <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(SS_STRATA, function (s) {
#     rs <- raw_slopes[[p]][[s]]
#     rs$q <- cut(rs$raw_slope, quantile(rs$raw_slope), include.lowest = T,
#                 labels = paste0("q", 1:4))
#     df <- long_adiposity[[p]]
#     df <- subset(df, df$eid %in% rs$eid)
#     df$q <- rs$q[match(df$eid, rs$eid)]
#     
#     # Calculate mean and SE in each 5-year interval within each quartile
#     age_bin_cuts <- seq(20, 80, by = 5)
#     df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
#     plot_df <- df %>% group_by(q, age_bin) %>% 
#       summarise(count = n(),
#                 mean_value = mean(value),
#                 se_value = sd(value)/sqrt(count))
#     
#     # Plot 
#     p <- ggplot(plot_df, aes(x = age_bin, y = mean_value,
#                              group = q, color = q, fill = q)) +
#       geom_point() +
#       geom_path() +
#       geom_ribbon(aes(ymin = mean_value - se_value, 
#                       ymax = mean_value + se_value),
#                   alpha = 0.2) +
#       scale_fill_brewer(palette = "Set1") +
#       scale_color_brewer(palette = "Set1") +
#       labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
#            title = s)
#     return (p)
#   })
#   names(res) <- SS_STRATA
#   pdf(paste0("plots/raw_slopes/mean_trajectories_quartile_", p, ".pdf"),
#       onefile = T)
#   print(res)
#   dev.off()
#   return (res)
# })
# names(rs_quartile_trajectories) <- PHENOTYPES
# 
# # Plot correlations between raw slope and model covariates ----
# 
# rs_correlations <- lapply(PHENOTYPES, function (p) {
#   
#   res <- lapply(SC_STRATA, function (s) {
#     df <- raw_slopes[[p]][[s]][, c("eid", "raw_slope", 
#                              "sex", "ancestry", 
#                              "FU_n", "FUyrs", 
#                              "baseline_age", "height", 
#                              "baseline_BMI", "baseline_trait")]
#     # Sample 1000 IDs to plot because otherwise too many to interpret
#     SAMPLE_IDS <- sample(df$eid, min(1000, nrow(df)), replace = F)
#     # But calculate correlations from all
#     plot_df <- pivot_longer(df, cols = c(FU_n, FUyrs, baseline_age, height,
#                                          baseline_BMI, baseline_trait),
#                             names_to = "covariate", values_to = "covariate_value")
#     point_plots <- subset(plot_df, plot_df$eid %in% SAMPLE_IDS)
#     
#     # Plot correlations 
#     res <- ggplot(plot_df, aes(x = covariate_value, y = raw_slope,
#                                color = sex, fill = sex)) +
#       facet_wrap(~covariate, ncol = 2, scales = "free_x") + 
#       geom_point(data = point_plots, 
#                  aes(x = covariate_value, y = raw_slope, color = sex),
#                  alpha = 0.5) +
#       geom_smooth(method = "lm") +
#       labs(x = "Covariate value", y = "Raw slope", 
#            title = paste(unique(df$ancestry), "ancestry"))
#     
#     return (res)
#   })
#   names(res) <- SC_STRATA
#   pdf(paste0("plots/raw_slopes/covariate_correlations_", p, ".pdf"),
#       onefile = T)
#   print(res)
#   dev.off()
#   return (res)
# })
# names(rs_correlations) <- PHENOTYPES
# 
# 
# # DEBUGGING: # Plot characteristics of q1 vs q4 by raw slope, white ----
# 
# rs_white_descriptive <- lapply(PHENOTYPES, function (p) {
#   df <- lapply(c("white_F", "white_M"), function (s) {
#     rs <- raw_slopes[[p]][[s]]
#     rs$q <- cut(rs$raw_slope, quantile(rs$raw_slope), include.lowest = T,
#                 labels = paste0("q", 1:4))
#     rs <- subset(rs, rs$q %in% c("q1", "q4"))
#     return (rs)
#   })
#   df <- bind_rows(df)
#   
#   # Number of follow-up measures
#   FU_n <- ggplot(df, aes(x = q, y = FU_n, fill = q)) +
#     facet_wrap(~sex, ncol = 2) +
#     geom_boxplot(position = position_dodge(1)) +
#     scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
#     labs(x = "raw slope quartile", y = "# repeat measures") +
#     theme(legend.position = "none")
#   
#   # Number of years of follow up
#   FUyrs <- ggplot(df, aes(x = q, y = FUyrs, fill = q)) +
#     facet_wrap(~sex, ncol = 2) +
#     geom_boxplot(position = position_dodge(1)) +
#     scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
#     labs(x = "raw slope quartile", y = "# follow-up years") +
#     theme(legend.position = "none")
#   
#   # Baseline age
#   bl_age <- ggplot(df, aes(x = q, y = baseline_age, fill = q)) +
#     facet_wrap(~sex, ncol = 2) +
#     geom_boxplot(position = position_dodge(1)) +
#     scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
#     labs(x = "raw slope quartile", y = "baseline age") +
#     theme(legend.position = "none")
#   
#   # Baseline BMI
#   bl_BMI <- ggplot(df, aes(x = q, y = baseline_BMI)) +
#     facet_wrap(~sex, ncol = 2) +
#     geom_violin(aes(fill = q), position = position_dodge(1)) +
#     geom_boxplot(width = 0.1) + 
#     scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
#     labs(x = "raw slope quartile", y = "baseline BMI (kg/m2)") +
#     theme(legend.position = "none")
#   
#   # Baseline adiposity trait
#   bl_trait <- ggplot(df, aes(x = q, y = baseline_trait)) +
#     facet_wrap(~sex, ncol = 2) +
#     geom_violin(aes(fill = q), position = position_dodge(1)) +
#     geom_boxplot(width = 0.1) + 
#     scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
#     labs(x = "raw slope quartile", y = "baseline trait") +
#     theme(legend.position = "none")
#   
#   pdf(paste0("plots/raw_slopes/white_q1_v_q4_covariate_distributions_", p, ".pdf"),
#       onefile = T)
#   print(FU_n)
#   print(FUyrs)
#   print(bl_age)
#   print(bl_BMI)
#   print(bl_trait)
#   dev.off()
#   return ()
# })
# 
# # Plot trajectories of sample individuals at the tails of adj slope distribution ----
# 
# as_tail_trajectories <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(STRATA, function (s) {
#     rs_df <- adj_slopes[[p]][[s]]
#     top_ids <- rs_df$eid[rs_df$residual > quantile(rs_df$residual, 0.95)]
#     top_ids <- sample(top_ids, min(length(top_ids), 5), replace = F)
#     bottom_ids <- rs_df$eid[rs_df$residual < quantile(rs_df$residual, 0.05)]
#     bottom_ids <- sample(bottom_ids, min(length(bottom_ids), 5), replace = F)
#     
#     plot_df <- long_adiposity[[p]]
#     plot_df <- subset(plot_df, 
#                       plot_df$eid %in% top_ids | plot_df$eid %in% bottom_ids)
#     plot_df$tail <- ifelse(plot_df$eid %in% top_ids, "top", "bottom")
#     
#     p <- ggplot(plot_df, aes(x = age_event, y = value, 
#                              group = eid, color = tail)) +
#       geom_line() +
#       geom_point() +
#       scale_color_manual(values = c("top" = "#984EA3", "bottom" = "#E41A1C"), 
#                          guide = F) +
#       labs(x = "age (years)", y = "Adiposity trait value", 
#            title = s)
#     return (p)
#   })
#   names(res) <- STRATA
#   pdf(paste0("plots/adjusted_slopes/tail_trajectories_", p, ".pdf"),
#       onefile = T)
#   print(res)
#   dev.off()
#   return (res)
# })
# names(as_tail_trajectories) <- PHENOTYPES
# 
# # Plot mean (binned) trajectories by adj slope quartile ----
# 
# as_quartile_trajectories <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(STRATA, function (s) {
#     rs <- adj_slopes[[p]][[s]]
#     rs$q <- cut(rs$residual, quantile(rs$residual), include.lowest = T,
#                 labels = paste0("q", 1:4))
#     df <- long_adiposity[[p]]
#     df <- subset(df, df$eid %in% rs$eid)
#     df$q <- rs$q[match(df$eid, rs$eid)]
#     
#     # Calculate mean and SE in each 5-year interval within each quartile
#     age_bin_cuts <- seq(20, 80, by = 5)
#     df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
#     plot_df <- df %>% group_by(q, age_bin) %>% 
#       summarise(count = n(),
#                 mean_value = mean(value),
#                 se_value = sd(value)/sqrt(count))
#     
#     # Plot 
#     p <- ggplot(plot_df, aes(x = age_bin, y = mean_value,
#                              group = q, color = q, fill = q)) +
#       geom_point() +
#       geom_path() +
#       geom_ribbon(aes(ymin = mean_value - se_value, 
#                       ymax = mean_value + se_value),
#                   alpha = 0.2) +
#       scale_fill_brewer(palette = "Set1") +
#       scale_color_brewer(palette = "Set1") +
#       labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
#            title = s)
#     return (p)
#   })
#   names(res) <- STRATA
#   pdf(paste0("plots/adjusted_slopes/mean_trajectories_quartile_", p, ".pdf"),
#       onefile = T)
#   print(res)
#   dev.off()
#   return (res)
# })
# names(as_quartile_trajectories) <- PHENOTYPES
# 
# # Plot mean (binned) trajectories by gainer status ----
# 
# gainer_status_trajectories <- lapply(PHENOTYPES, function (p) {
#   res <- lapply(STRATA, function (s) {
#     rs <- adj_slopes[[p]][[s]]
#     df <- long_adiposity[[p]]
#     df <- subset(df, df$eid %in% rs$eid)
#     df$gainer <- rs$gainer[match(df$eid, rs$eid)]
#     
#     # Calculate mean and SE in each 5-year interval within gainers/non-gainers
#     age_bin_cuts <- seq(20, 80, by = 5)
#     df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
#     plot_df <- df %>% group_by(gainer, age_bin) %>% 
#       summarise(count = n(),
#                 mean_value = mean(value),
#                 se_value = sd(value)/sqrt(count))
#     
#     plot_df$ancestry <- strsplit(s, "_")[[1]][[1]]
#     plot_df$sexplot <- strsplit(s, "_")[[1]][[2]]
#     return (plot_df)
#   })
#   res <- bind_rows(res)
#   res <- split(res, res$ancestry)
#   res_plots <- lapply(res, function (a) {
#     # Plot 
#     res_plots <- ggplot(a, aes(x = age_bin, y = mean_value,
#                              group = gainer, color = gainer, fill = gainer)) +
#       facet_wrap(~sexplot, nrow = 3) +
#       geom_point() +
#       geom_path() +
#       geom_ribbon(aes(ymin = mean_value - se_value, 
#                       ymax = mean_value + se_value),
#                   alpha = 0.2) +
#       scale_color_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C"), 
#                          guide = F) +
#       scale_fill_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C"), 
#                          guide = F) +
#       labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
#            title = unique(a$ancestry))
#     return (res_plots)
#   })
#   pdf(paste0("plots/adjusted_slopes/mean_trajectories_gainer_status_", p, ".pdf"),
#       onefile = T)
#   print(res_plots)
#   dev.off()
#   return (res_plots)
# })
# names(gainer_status_trajectories) <- PHENOTYPES
# 
