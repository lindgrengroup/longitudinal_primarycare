# Author: Samvida S. Venkatesh
# Date: 02/03/21

library(tidyverse)
#library(treemap)
library(RColorBrewer)
theme_set(theme_bw())

# Read stratified raw slope, adiposity, and covariate data ----

raw_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_raw_slopes.rds")
adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/QCd_adiposity.rds")
PHENOTYPES <- names(raw_slopes)
STRATA <- unique(unlist(lapply(raw_slopes, function (x) names(x) )))

# Construct demographic characteristics tables -----

demo_table <- lapply(PHENOTYPES, function (p) {
  pheno <- lapply(raw_slopes[[p]], function (s) {
    s[, c("eid", "sex", "ancestry", "FU_n", "FUyrs", "baseline_age", "height",
          "baseline_BMI")]
  })
  df <- bind_rows(pheno)
  summary_df <- df %>% group_by(sex, ancestry) %>% 
    summarise(count = n(), 
              mean_FU_n = mean(FU_n), se_FU_n = sd(FU_n)/sqrt(count),
              mean_FUyrs = mean(FUyrs), se_FUyrs = sd(FUyrs)/sqrt(count),
              mean_bl_age = mean(baseline_age), 
              se_bl_age = sd(baseline_age)/sqrt(count),
              median_height = median(height), 
              iqr_height = paste(quantile(height, 0.25), quantile(height, 0.75),
                                 sep = ", "),
              median_bl_BMI = median(baseline_BMI), 
              iqr_bl_BMI = paste(quantile(baseline_BMI, 0.25), 
                                 quantile(baseline_BMI, 0.75),
                                 sep = ", ")) %>%
    ungroup() %>% mutate(percent = count / sum(count))
  
  adipo <- lapply(STRATA, function (s) {
    res <- adiposity[[p]][[s]][, c("eid", "age_event", "value")]
    res <- res %>% group_by(eid) %>% arrange(age_event) %>% 
      summarise(bl_value = first(value)) 
    res$ancestry <- strsplit(s, "_")[[1]][1]
    res$sex <- strsplit(s, "_")[[1]][2]
    return (res)
  })
  adipo <- bind_rows(adipo)
  adipo_summary <- adipo %>% group_by(sex, ancestry) %>% 
    summarise(median_bl_value = median(bl_value), 
              iqr_bl_value = paste(quantile(bl_value, 0.25), 
                                 quantile(bl_value, 0.75),
                                 sep = ", "))
  summary_df <- merge(summary_df, adipo_summary, by = c("sex", "ancestry"))
  
  write.table(summary_df, paste0("results/descriptive_factors_", p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  return (summary_df)
})

# Plot: ancestry and sex population counts ----

# demo_table <- lapply(PHENOTYPES, function (p) {
#   res <- read.table(paste0("results/descriptive_factors/descriptive_factors_", p, ".txt"),
#                     sep = "\t", header = T, stringsAsFactors = F)
# })
# names(demo_table) <- PHENOTYPES

popn_ancestry_treemap <- lapply(PHENOTYPES, function (p) {
  df <- demo_table[[p]] %>% group_by(ancestry) %>% summarise(size = sum(count))
  tiff(filename = paste0("plots/descriptive_factors/popn_ancestry_", p, ".tiff"),
       height = 10, width = 10, units = "cm",
       res = 600)
  treemap(df, 
          index = "ancestry", vSize = "size", type = "index",
          palette = "Dark2", border.col = c("black"), border.lwds = 1,
          title = "", fontsize.labels = 0)
  dev.off()
  return ()
})

popn_sex_pies <- lapply(PHENOTYPES, function (p) {
  df <- demo_table[[p]]
  for (a in unique(df$ancestry)) {
    res <- ggplot(subset(df, df$ancestry == a), 
                  aes(x = "", y = count, fill = sex)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0) +
      scale_color_brewer(palette = "Dark2") +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.title.y = element_blank(), panel.border = element_blank(), 
            panel.grid = element_blank(), axis.ticks = element_blank(),
            legend.position = "none")
    tiff(filename = paste0("plots/descriptive_factors/popn_sex_", p, "_", 
                           a, ".tiff"),
         height = 10, width = 10, units = "cm",
         res = 600)
    print(res)
    dev.off()
  }
  return ()
})

# Box and violin plots of distributions ----

plot_distributions <- lapply(PHENOTYPES, function (p) {
  pheno <- lapply(raw_slopes[[p]], function (s) {
    s[, c("eid", "sex", "ancestry", "FU_n", "FUyrs", "baseline_age", "height",
          "baseline_BMI")]
  })
  df <- bind_rows(pheno)
  
  # Number of follow-up measures
  FU_n <- ggplot(df, aes(x = ancestry, y = FU_n, fill = ancestry)) +
    facet_wrap(~sex, nrow = 2) +
    geom_boxplot(position = position_dodge(1)) +
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "# repeat measures") +
    theme(legend.position = "none")
  
  # Number of years of follow up
  FUyrs <- ggplot(df, aes(x = ancestry, y = FUyrs, fill = ancestry)) +
    facet_wrap(~sex, nrow = 2) +
    geom_boxplot(position = position_dodge(1)) +
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "follow-up length (years)") +
    theme(legend.position = "none")
  
  # Baseline age
  bl_age <- ggplot(df, aes(x = ancestry, y = baseline_age)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(20, 80), breaks = seq(20, 80, by = 5)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "baseline age (years)") +
    theme(legend.position = "none")
  
  # Adult height
  med_height <- ggplot(df, aes(x = ancestry, y = height)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(110, 210), breaks = seq(110, 210, by = 10)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "median adult height (cm)") +
    theme(legend.position = "none")
  
  # Baseline BMI
  bl_BMI <- ggplot(df, aes(x = ancestry, y = baseline_BMI)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(10, 65), breaks = seq(10, 65, by = 5)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "baseline BMI (kg/m2)") +
    theme(legend.position = "none")
  
  # Baseline adiposity trait
  adipo <- lapply(STRATA, function (s) {
    res <- adiposity[[p]][[s]][, c("eid", "age_event", "value")]
    res <- res %>% group_by(eid) %>% arrange(age_event) %>% 
      summarise(bl_value = first(value)) 
    res$ancestry <- strsplit(s, "_")[[1]][1]
    res$sex <- strsplit(s, "_")[[1]][2]
    return (res)
  })
  adipo <- bind_rows(adipo)
  bl_adipo <- ggplot(adipo, aes(x = ancestry, y = bl_value)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "baseline adiposity trait") +
    theme(legend.position = "none")
  
  pdf(paste0("plots/descriptive_factors/covariate_distributions_", p, ".pdf"),
      onefile = T)
  print(FU_n)
  print(FUyrs)
  print(bl_age)
  print(med_height)
  print(bl_BMI)
  print(bl_adipo)
  dev.off()
  return ()
})

# Summary plot of baseline age, FU length, FU measures -----

# demo_table <- lapply(PHENOTYPES, function (p) {
#   res <- read.table(paste0("results/descriptive_factors/descriptive_factors_", 
#                            p, ".txt"),
#                     sep = "\t", header = T, stringsAsFactors = F)
# })
# names(demo_table) <- PHENOTYPES

demo_table <- lapply(PHENOTYPES, function (p) {
  res <- demo_table[[p]]
  res$adipo_trait <- p
  return (res)
})
demo_table <- bind_rows(demo_table)
# min-max normalise line thickness from 1:5
demo_table$line_thick <- 1 + 
  4*(demo_table$mean_FU_n - min(demo_table$mean_FU_n)) / 
  (max(demo_table$mean_FU_n) - min(demo_table$mean_FU_n))

pdf("plots/descriptive_factors/summary_FU.pdf")
ggplot(demo_table) +
  facet_wrap(~ancestry+sex, ncol = 2) +
  geom_segment(aes(x = mean_bl_age, xend = mean_bl_age + mean_FUyrs,
                   y = adipo_trait, yend = adipo_trait, 
                   size = mean_FU_n,
                   color = ancestry)) +
  scale_y_discrete(limits = c("WHR", "waist", "weight", "BMI"),
                   labels = c("WHR", "WC", "weight", "BMI")) +
  scale_x_continuous(limits = c(40, 65), breaks = seq(40, 65, by = 5)) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Age (years)", y = "Adiposity trait") +
  theme(legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank())
dev.off()

# Plot binned trajectories ----

# Calculate mean and SE in each 5-year interval in each strata

age_bin_cuts <- seq(20, 80, by = 5)

binned_adipo <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    df <- adiposity[[p]][[s]]
    df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
    res <- df %>% group_by(age_bin) %>% summarise(count = n(),
                                                  mean_value = mean(value),
                                                  se_value = sd(value)/sqrt(count))
    res$ancestry <- strsplit(s, "_")[[1]][1]
    res$sex <- strsplit(s, "_")[[1]][2]
    return (res)
  })
  res <- bind_rows(res)
  return (res)
})
names(binned_adipo) <- PHENOTYPES

traj_plots <- lapply(PHENOTYPES, function (p) {
  res <- ggplot(binned_adipo[[p]], aes(x = age_bin, y = mean_value,
                                group = ancestry,
                                color = ancestry, fill = ancestry)) +
    facet_wrap(~sex, nrow = 2) +
    geom_point() +
    geom_path() +
    geom_ribbon(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                alpha = 0.2) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
         title = p)
  return (res)
})
names(traj_plots) <- PHENOTYPES

pdf("plots/trajectories/binned_raw_data.pdf", onefile = T)
print(traj_plots)
dev.off()
# Construct raw slope summary table ----

rs_table <- lapply(PHENOTYPES, function (p) {
  pheno <- lapply(raw_slopes[[p]], function (s) {
    s[, c("eid", "sex", "ancestry", "raw_slope")]
  })
  df <- bind_rows(pheno)
  summary_df <- df %>% group_by(sex, ancestry) %>% 
    summarise(median_slope = median(raw_slope), 
              iqr_slope = paste(quantile(raw_slope, 0.25), 
                                quantile(raw_slope, 0.75),
                                 sep = ", "))
  write.table(summary_df, paste0("results/raw_slope_summaries_", p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  return (summary_df)
})

# Plot violin distribution of raw slopes ----

rs_distributions <- lapply(PHENOTYPES, function (p) {
  pheno <- lapply(raw_slopes[[p]], function (s) {
    s[, c("eid", "sex", "ancestry", "raw_slope")]
  })
  df <- bind_rows(pheno)
  res <- ggplot(df, aes(x = ancestry, y = raw_slope)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "raw slope", title = p) +
    theme(legend.position = "none")
  return (res)
})
pdf(paste0("plots/raw_slopes/raw_slopes_distributions.pdf"),
    onefile = T)
print(rs_distributions)
dev.off()
