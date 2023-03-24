# Author: Samvida S. Venkatesh
# Date: 23/03/23

library(tidyverse)
theme_set(theme_bw())

# Read data ----

maindat_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data" # REDACTED
gpdat_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity"
gen_resources_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources" # REDACTED
outfile_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_rs429358" # REDACTED

# Genotypes / dosages at rs429358
apoe_dosages <- read.table(paste0(gpdat_path, "/sample_variant_counts/rs429358_dosages.txt"),
                           sep = " ", header = T, stringsAsFactors = F)
# Remove first row, which contains info on type of column and columns 
# 2, 3, 4 (ID repeat, missingness, sex)
apoe_dosages <- apoe_dosages[-1, c(1, 5)]
colnames(apoe_dosages) <- c("eid", "dosage")
apoe_dosages$eid <- as.character(apoe_dosages$eid)

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

covars <- readRDS(paste0(maindat_path, "/covariates.rds"))[PHENOTYPES]
general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Wrangle data to test associations ----

full_dat <- lapply(covars, function (df) {
  df$eid <- as.character(df$eid)
  res <- df %>% inner_join(apoe_dosages,
                           by = "eid")
  res$sex <- general_covars$sex[match(res$eid, general_covars$eid)]
  res$dosage <- as.numeric(res$dosage)
  return (res)
})
names(full_dat) <- PHENOTYPES

# Test associations ----

res_table <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    sub_df <- full_dat[[p]]
    if (sx != "sex_comb") {
      sub_df <- sub_df %>% filter(sex == sx)
    }
    
    # SNP association with length of follow-up (years)
    mod_fuyrs <- lm(FUyrs ~ dosage, data = sub_df)
    fuyrs_summ <- as.data.frame(t(summary(mod_fuyrs)$coefficients["dosage", c(1,2,4)]))
    fuyrs_summ$outcome <- "FU_yrs"
    # SNP association with number of follow-up measures
    mod_fun <- lm(FU_n ~ dosage, data = sub_df)
    fun_summ <- as.data.frame(t(summary(mod_fun)$coefficients["dosage", c(1,2,4)]))
    fun_summ$outcome <- "FU_n"
    
    res <- bind_rows(fuyrs_summ, fun_summ) %>%
      mutate(sex_strata = sx)
    
    return (res)
  })
  res_list <- bind_rows(res_list) %>%
    mutate(pheno = p)
  return (res_list)
})
res_table <- bind_rows(res_table)

write.table(res_table,
            paste0(outfile_path, "/assocn_with_followup.txt"),
            sep = "\t", quote = F, row.names = F)

# Plot associations ----

custom_two_diverge <- c("#D35C79","#009593")
names(custom_two_diverge) <- c("F", "M")

lapply(PHENOTYPES, function (p) {
  plot_dat <- full_dat[[p]] %>%
    filter(sex %in% c("F", "M"))
  plot_dat$dosage <- factor(plot_dat$dosage, levels = c("0", "1", "2"))
  
  fuyrs_plot <- ggplot(plot_dat,
                       aes(x = dosage, y = FUyrs, 
                           color = sex, fill = sex)) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(aes(group = interaction(sex, dosage)),
                 width = 0.1, position = position_dodge(width = 0.9),
                 colour = "black", fill = "white") +
    scale_fill_manual(values = custom_two_diverge, guide = "none") +
    scale_color_manual(values = custom_two_diverge, guide = "none") +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 8))
  ggsave(paste0(outfile_path, "/plot_fuyrs_", p, ".png"),
         fuyrs_plot, units = "in", height = 3, width = 6.5)
  
  fun_plot <- ggplot(plot_dat,
                     aes(x = dosage, y = FU_n, 
                         color = sex, fill = sex)) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(aes(group = interaction(sex, dosage)),
                 width = 0.1, position = position_dodge(width = 0.9),
                 colour = "black", fill = "white") +
    scale_fill_manual(values = custom_two_diverge, guide = "none") +
    scale_color_manual(values = custom_two_diverge, guide = "none") +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 8))
  ggsave(paste0(outfile_path, "/plot_fu_n_", p, ".png"),
         fun_plot, units = "in", height = 3, width = 6.5)
  
  fun_cut_plot <- ggplot(plot_dat,
                     aes(x = dosage, y = FU_n, 
                         color = sex, fill = sex)) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(aes(group = interaction(sex, dosage)),
                 width = 0.1, position = position_dodge(width = 0.9),
                 colour = "black", fill = "white") +
    scale_fill_manual(values = custom_two_diverge, guide = "none") +
    scale_color_manual(values = custom_two_diverge, guide = "none") +
    scale_y_continuous(limits = c(0, 50)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_text(size = 8))
  ggsave(paste0(outfile_path, "/plot_fu_n_cut_", p, ".png"),
         fun_cut_plot, units = "in", height = 3, width = 6.5)
})


