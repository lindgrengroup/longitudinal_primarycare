# Author: Samvida S. Venkatesh
# Date: 25/11/2020

# Remove extreme BMI values based on full UKBIOBANK data
# Values +/- 10% more extreme than the min and max UKBIOBANK BMI are removed as noise

library(tidyverse)
library(dplyr)
library(RColorBrewer)
theme_set(theme_bw())
library(plotly)

# Calculate thresholds ----

pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

ukbb_bmi <- pheno[, c("f.eid", 
                      "f.21001.0.0", "f.21001.1.0", "f.21001.2.0")]

colnames(ukbb_bmi) <- c("eid", "BMI_1", "BMI_2", "BMI_3")
ukbb_bmi <- pivot_longer(ukbb_bmi, starts_with("BMI"),
                         names_to = "Event",
                         values_to = "BMI",
                         values_drop_na = T)

sink("/well/lindgren/UKBIOBANK/samvida/BMI/UKBIOBANK_BMI_summary.txt")
print(summary(ukbb_bmi$BMI))
print(paste("MIN_THRESHOLD:", 0.9*min(ukbb_bmi$BMI, na.rm = T)))
print(paste("MAX_THRESHOLD:", 1.1*max(ukbb_bmi$BMI, na.rm = T)))
sink()

# Read data ----

pc_data <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/bmi_primary_care_annotated.txt",
                      sep = "\t", header = T)

# Look at BMI trajectories for individuals in the top and bottom 0.1% of UKBIOBANK
# to evaluate cutoffs (proposed: 10.91 - 82.15)

UKBB_BMIs <- distinct(pc_data, eid, sex, mean_UKBB_BMI)
qs <- UKBB_BMIs %>% group_by(sex) %>% 
  mutate(q_low = quantile(mean_UKBB_BMI, 0.001, na.rm = T)[[1]],
            q_high = quantile(mean_UKBB_BMI, 0.999, na.rm = T)[[1]]) %>%
  mutate(light = mean_UKBB_BMI < q_low,
         heavy = mean_UKBB_BMI > q_high)

light <- pc_data[pc_data$eid %in% qs$eid[which(qs$light)], ]
heavy <- pc_data[pc_data$eid %in% qs$eid[which(qs$heavy)], ]

tmp <- ggplot(light, aes(x = age_years, y = primarycare_BMI, col = sex)) +
  geom_point() +
  geom_line() +
  labs(x = "Age (years)", y = "BMI")

htmlwidgets::saveWidget(as_widget(ggplotly(tmp)), 
                        "/well/lindgren/UKBIOBANK/samvida/BMI/plots/trajectories_light.html")


tmp <- ggplot(heavy, aes(x = age_years, y = primarycare_BMI, col = sex)) +
  geom_point() +
  geom_line() +
  labs(x = "Age (years)", y = "BMI")

htmlwidgets::saveWidget(as_widget(ggplotly(tmp)), 
                        "/well/lindgren/UKBIOBANK/samvida/BMI/plots/trajectories_heavy.html")

# Visualise thresholds ----

MIN <- 10.00
MAX <- 82.15

# Plot to visualise cut-off thresholds
pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/popn_outlier_thresholds.pdf")

ggplot(pc_data, aes(x = primarycare_BMI, color = sex, fill = sex)) +
  geom_density(alpha = 0.25, na.rm = T) +
  geom_vline(xintercept = MIN) +
  geom_vline(xintercept = MAX) +
  geom_vline(xintercept = 12.12, linetype = "dashed") + 
  geom_vline(xintercept = 74.68, linetype = "dashed") + 
  labs(x = "BMI", y = "Density", 
       title = "Primary Care BMI distribution and cutoff thresholds") +
  xlim(c(10, 100))

# Remove population-level outliers based on thresholds ----

pc_data$remove <- pc_data$primarycare_BMI < MIN | 
  pc_data$primarycare_BMI > MAX

# What do we lose? ----

removed <- pc_data[pc_data$remove, ]
cleaned <- pc_data[!pc_data$remove, ]

# Look at distribution of BMI observations that were removed for being too low
ggplot(subset(removed, removed$primarycare_BMI < MIN), 
       aes(x = primarycare_BMI, color = sex, fill = sex)) +
  geom_density(alpha = 0.25, na.rm = T) +
  xlim(0, MIN) +
  labs(x = "BMI", y = "Density", 
       title = "Distribution of BMI measures removed for being too low") 

# or too high
ggplot(subset(removed, removed$primarycare_BMI > MAX), 
       aes(x = primarycare_BMI, color = sex, fill = sex)) +
  geom_density(alpha = 0.25, na.rm = T) +
  xlim(MAX, max(removed$primarycare_BMI)) +
  labs(x = "BMI", y = "Density", 
       title = "Distribution of BMI measures removed for being too high") 

# Look at mean UKBIOBANK BMI of all individuals vs those we lose due to BMI measurement
# errors

all_indivs <- distinct(pc_data, eid, sex, mean_UKBB_BMI)
all_indivs$group <- "all"
removed_indivs <- all_indivs[!all_indivs$eid %in% cleaned$eid, ]
removed_indivs$group <- "removed"

table(removed_indivs$sex)

p <- bind_rows(all_indivs, removed_indivs)

ggplot(p, aes(x = mean_UKBB_BMI, col = group, fill = group)) +
  facet_wrap(~sex, nrow = 2) +
  geom_density(alpha = 0.25, na.rm = T) +
  xlim(10, 100) +
  labs(x = "BMI", y = "Density", 
       title = "Mean UKBIOBANK BMI distribution in all individuals vs those for whom
       no data remains after error removal") 

dev.off()

write.table(cleaned[, -8], 
            "/well/lindgren/UKBIOBANK/samvida/BMI/bmi_primary_care_errors_removed.txt",
            sep = "\t", quote = F, row.names = F)