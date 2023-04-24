# Author: Samvida S. Venkatesh
# Date: 25/03/23

library(tidyverse)
theme_set(theme_bw())

# Read files ----

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
outpath <- "" # REDACTED

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

dat <- readRDS(paste0(mainpath, "/indiv_qcd_data.rds"))[PHENOTYPES]
covars <- readRDS(paste0(mainpath, "/covariates.rds"))[PHENOTYPES]

general_covars <- read.table(paste0(gen_resources_path, "/220504_QCd_demographic_covariates.txt"),
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

# Wrangle data ----

full_dat <- lapply(PHENOTYPES, function (p) {
  df <- dat[[p]]
  df$eid <- as.character(df$eid)
  df$sex <- general_covars$sex[match(df$eid, general_covars$eid)]
  
  covars$eid <- as.character(covars$eid)
  df <- df %>% left_join(covars[[p]], by = "eid")
  return (df)
})
names(full_dat) <- PHENOTYPES

# Characterise variance within each strata ----

custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

plotVar <- function (vardf, pop_var, sx) {
  ggplot(vardf,
         aes(x = var_value)) +
    geom_density(fill = custom_three_diverge[sx],
                 colour = custom_three_diverge[sx],
                 alpha = 0.7) +
    geom_vline(xintercept = pop_var, linetype = "dashed") +
    scale_x_continuous(trans = "log10",
                       guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))
}

var_calc <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    cat(paste0("Pheno: ", p, ", sex strata: ", sx), "\n")
    df <- full_dat[[p]]
    if (sx != "sex_comb") df <- df %>% filter(sex == sx)
    
    # Individual-level
    indiv_var <- df %>%
      group_by(eid) %>%
      summarise(var_value = var(value, na.rm = T))
    cat(paste0("Mean intra-individual variance: ", mean(indiv_var$var_value,
                                                        na.rm = T)), "\n")
    
    # Population-level
    pop_var <- var(df$value, na.rm = T)
    cat(paste0("Inter-individual variance: ", pop_var), "\n")
    
    # Save plot
    ggsave(paste0(outpath, "/plots/intra_indiv_variance_", p, "_", sx,
                  ".png"),
           plotVar(indiv_var, pop_var, sx),
           units = "in", height = 2, width = 2)
    return (indiv_var)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(var_calc) <- PHENOTYPES

# Correlation between variance and follow-up ----

plotCorrs <- function (df, fu_var) {
  ggplot(df, aes(x = var_value,
                 y = !!as.symbol(fu_var), colour = sex, fill = sex)) +
    geom_point(size = 0.3, alpha = 0.3) +
    scale_x_continuous(trans = "log10",
                       guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))
}

corrvar_FU <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    cat(paste0("Pheno: ", p, ", sex strata: ", sx), "\n")
    
    df <- var_calc[[p]][[sx]]
    df$eid <- as.character(df$eid)
    covars[[p]]$eid <- as.character(covars[[p]]$eid)
    
    df <- left_join(df, covars[[p]],
                    by = "eid")
    df$sex <- general_covars$sex[match(df$eid, general_covars$eid)]
    df <- df[complete.cases(df), ]
    cat(paste0("Variance correlation with FUyrs: ", cor(df$var_value,
                                                        df$FUyrs)), "\n")
    cat(paste0("Variance correlation with FU_n: ", cor(df$var_value,
                                                        df$FU_n)), "\n")
    
    if (sx == "sex_comb") {
      ggsave(paste0(outpath, "/plots/variance_vs_FUyrs_", p, "_", sx,
                    ".png"),
             plotCorrs(df, "FUyrs"),
             units = "in", height = 2, width = 2)
      ggsave(paste0(outpath, "/plots/variance_vs_FU_n_", p, "_", sx,
                    ".png"),
             plotCorrs(df, "FU_n"),
             units = "in", height = 2, width = 2)
    }
    
  })
})

