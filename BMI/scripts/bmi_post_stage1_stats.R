# Author: Samvida S. Venkatesh
# Date: 26/11/20

library(tidyverse)
library(dplyr)
theme_set(theme_bw())
library(RColorBrewer)

BMI <- read.table("/well/lindgren/UKBIOBANK/samvida/BMI/bmi_primary_care_errors_removed.txt",
                  sep = "\t", header = T)

# Data summaries ----

table(BMI$sex)

indivs <- BMI %>% group_by(eid, sex) %>% count()
table(indivs$sex)

p <- data.frame(type = c("observations", "observations",
                         "individuals", "individuals"),
                sex = c("F", "M", "F", "M"),
                counts = c(785539, 645940, 114434, 94430))

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/post_stage1/summaries.pdf")

# Number of individuals and observations
ggplot(p, aes(x = type, y = counts, color = sex, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "", y = "Count")

# Distribution of BMI values
ggplot(BMI, aes(x = primarycare_BMI, color = sex, fill = sex)) +
  geom_density(alpha = 0.25) +
  labs(x = "BMI", y = "Density")

# Distribution of BMI by sex (violins)
ggplot(BMI, aes(x = sex, y = primarycare_BMI)) +
  geom_jitter(size = 0.1, color = "grey") +
  geom_violin(aes(fill = sex), trim = F) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  labs(x = "Sex", y = "BMI")

# Distribution of BMI by age (boxes)
breakpoints <- c(18, seq(30, 80, by = 10))
names <- c("(18-30]", "(30-40]", "(40-50]", "(50-60]", "(60-70]",
           "(70-80]")
BMI$age_10_bin <- cut(BMI$age_years, breaks = breakpoints, labels = names)

ggplot(BMI, aes(x = age_10_bin, y = primarycare_BMI, fill = sex)) +
  geom_boxplot() +
  labs(x = "Age (years)", y = "BMI") 

dev.off()

## Longitudinal events ----

indivs <- BMI %>% group_by(eid, sex, mean_UKBB_BMI) %>% 
  summarise(mean_primarycare_BMI = mean(primarycare_BMI), 
            n_obs = n())

# Group individuals by number of BMI measures (1, 2-5, 6-10, 11-15, 16+)
# Create categories for number of measurements
breakpoints <- c(-Inf, 1, 5, 10, 15, Inf)
names <- c("1", "2-5", "6-10", "11-15", "16+")
indivs$nmeasures_bin <- cut(indivs$n_obs, breaks = breakpoints, labels = names)

# Create categories for BMI (underweight, normal, overweight, etc.)
# based on WHO guidelines - how many individuals are in each category?

breakpoints <- c(-Inf, 18.5, 25, 30, 35, 40, Inf)
names <- c("Underweight (< 18.5)", "Normal weight [18.5 - 25)", 
           "Pre-obesity [25 - 30)", "Obesity Class I [30 - 35)", 
           "Obesity Class II [35 - 40)", "Obesity Class III (>= 40)")
indivs$BMI_primarycare_class <- cut(indivs$mean_primarycare_BMI, 
                                 breaks = breakpoints, 
                                 include.lowest = T, right = F,
                                 labels = names)

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/post_stage1/longitudinal_summaries.pdf")

# Histogram of number of BMI measures
ggplot(indivs, aes(x = n_obs, fill = sex, colour = sex)) +
  geom_histogram(alpha = 0.25, position = "identity") +
  labs(x = "Number of BMI measures", y = "Number of individuals")
# zoom in to where the most data is
ggplot(indivs, aes(x = n_obs, fill = sex, colour = sex)) +
  geom_histogram(alpha = 0.25, position = "identity") +
  xlim(c(0, 30)) +
  labs(x = "Number of BMI measures", y = "Number of individuals")

# BMI distribution stratified by # times BMI is measured
med_grouped <- indivs %>% group_by(nmeasures_bin) %>% 
  summarise(med_in_group = median(mean_primarycare_BMI))

ggplot(indivs, aes(x = mean_primarycare_BMI, 
                   fill = nmeasures_bin, color = nmeasures_bin)) +
  geom_density(alpha = 0.25) +
  geom_vline(data = med_grouped, aes(xintercept = med_in_group, 
                                      color = nmeasures_bin)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(color = "# BMI measures", fill = "# BMI measures", 
       x = "Mean BMI")

# Number of times BMI is measured stratified by obesity status
ggplot(indivs, aes(x = BMI_primarycare_class, y = n_obs, fill = sex)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "BMI Class", 
       y = "# BMI measures")

## Age - number of measures - BMI relationship ----

# Group measurement age into 10-year bins
breakpoints <- c(18, seq(30, 80, by = 10))
names <- c("(18-30]", "(30-40]", "(40-50]", "(50-60]", "(60-70]",
           "(70-80]")
BMI$age_10_bin <- cut(BMI$age_years, breaks = breakpoints, labels = names)

test_age_effect <- group_split(BMI, age_10_bin)

# Create categories for number of measurements vs mean BMI within each age cut
breakpoints <- c(-Inf, 1, 5, 10, 15, Inf)
names <- c("1", "2-5", "6-10", "11-15", "16+")
test_age_effect <- lapply(test_age_effect, function (df) {
  df <- df %>% group_by(eid, sex, age_10_bin, mean_UKBB_BMI) %>% 
    summarise(mean_primarycare_BMI = mean(primarycare_BMI), 
              n_obs = n())
  df$nmeasures_bin <- cut(df$n_obs, breaks = breakpoints, labels = names)
  return (df)
})

test_age_effect <- bind_rows(test_age_effect)

ggplot(test_age_effect, aes(x = mean_primarycare_BMI, 
                   fill = nmeasures_bin, color = nmeasures_bin)) +
  facet_wrap(~age_10_bin, ncol = 2) +
  geom_density(alpha = 0.25) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(color = "# BMI measures", fill = "# BMI measures", 
       x = "Mean BMI")

dev.off()

