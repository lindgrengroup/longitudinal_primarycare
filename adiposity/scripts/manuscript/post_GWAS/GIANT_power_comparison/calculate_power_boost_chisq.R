# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(argparse)
library(tidyverse)
theme_set(theme_bw())

# Read files ----

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Which strata we are assessing")
parser$add_argument("--sample_size_gp", required=TRUE,
                    help = "Sample size in GP data")
args <- parser$parse_args()

STRATA <- args$strata

postgwas_path <- "" # REDACTED
sumstats_path <- "" # REDACTED
giant_sumstats_path <- "" # REDACTED

power_log <- paste0(postgwas_path, "/", STRATA, "/power_boost_log.txt")
power_plot <- paste0(postgwas_path, "/", STRATA, "/power_boost_plot.tiff")
power_res <- paste0(postgwas_path, "/", STRATA, "/power_boost_snp_comparisons.txt")

gp_dat <- read.table(paste0(sumstats_path, "/", STRATA, "/tmp_gp_gwas.txt"), 
                     sep = "\t", header = T, stringsAsFactors = F)
giant_dat <- read.table(paste0(giant_sumstats_path, "/", STRATA, "/tmp_giant_meta.txt"), 
                     sep = " ", header = T, stringsAsFactors = F)

# Calculate chi-square statistics and merge files -----

gp_chisq <- gp_dat %>% distinct(SNP, .keep_all = T) %>%
  mutate(CHI_SQUARE_GP = (BETA/SE)^2,
         N_GP = as.numeric(args$sample_size_gp)) %>%
  rename(BETA_GP = BETA, SE_GP = SE) %>%
  select(SNP, BETA_GP, SE_GP, CHI_SQUARE_GP, N_GP)

giant_chisq <- giant_dat %>% distinct(SNP, .keep_all = T) %>%
  mutate(CHI_SQUARE_GIANT = (BETA/SE)^2) %>%
  rename(BETA_GIANT = BETA, SE_GIANT = SE, N_GIANT = N) %>%
  select(SNP, BETA_GIANT, SE_GIANT, CHI_SQUARE_GIANT, N_GIANT)

merged_dat <- full_join(gp_chisq, giant_chisq, by = "SNP") %>%
  mutate(CHI_SQ_RATIO = CHI_SQUARE_GP/CHI_SQUARE_GIANT,
         EXP_N_RATIO = N_GP/N_GIANT) %>%
  filter(!is.na(CHI_SQ_RATIO))

write.table(merged_dat, power_res, 
            sep = "\t", row.names = F, quote = F)

obs_med <- median(merged_dat$CHI_SQ_RATIO, na.rm = T)
exp_med <- median(merged_dat$EXP_N_RATIO, na.rm = T)
power_boost <- (obs_med / exp_med) * 100

# Plot results ----

res_plot <- ggplot(merged_dat, aes(x = CHI_SQUARE_GIANT, 
                                   y = CHI_SQUARE_GP)) +
  geom_point(size = 0.3) +
  geom_abline(intercept = 0, slope = exp_med, 
              colour = "blueviolet", linetype = "dashed") + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 8))

tiff(filename = power_plot, height = 5, width = 5, units = "cm",
     res = 300)
print(res_plot)
dev.off()  

# Print results to power boost log ----

sink(power_log, append = T)
cat(paste0("Across all significant SNPs, observed chi-square ratio: ", 
           obs_med, "\n"))
cat(paste0("Across all significant SNPs, expected chi-square ratio: ", 
           exp_med, "\n"))
cat(paste0("Power boost: ", 
           power_boost, "%", "\n"))
sink()
