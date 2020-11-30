# Author: Samvida S. Venkatesh
# Date: 21/08/20

library(tidyr)
library(dplyr)
library(eeptools)
library(ggplot2)
theme_set(theme_bw())
library(RColorBrewer)

# Read data ----

# Read gp_clinical file
gp_clinical <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt",
                          sep = "\t", header = T, comment.char = "$",
                          stringsAsFactors = F)

# Remove individuals who have withdrawn consent
withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20200820.csv", 
                        header = F)
gp_clinical <- subset(gp_clinical, !gp_clinical$eid %in% withdrawn$V1)

# Add information on age and sex from main UKBB phenotypes file
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

gp_clinical <- merge(gp_clinical, pheno[, c("f.eid", "f.31.0.0", 
                                            "f.34.0.0", "f.52.0.0")], 
                     by.x = "eid", by.y = "f.eid", all.x = T)
colnames(gp_clinical)[9:11] <- c("sex", "year_of_birth", "month_of_birth")

gp_clinical$sex <- ifelse(gp_clinical$sex == 0, "F", "M")

# Assign date of birth (first of the month) based on month and year of birth 
gp_clinical$dob <- as.Date(paste(gp_clinical$year_of_birth, 
                                 gp_clinical$month_of_birth,
                                 1, sep = "-"))

# Extract read codes of interest ----

# Get relevant read codes

plasma_LH <- subset(gp_clinical, gp_clinical$read_3 == "XaELa" | 
                      gp_clinical$read_2 == "443f.")
plasma_LH$hormone <- "plasma LH"
serum_LH <- subset(gp_clinical, gp_clinical$read_3 %in% 
                     c("XM0lv", "X80Fu", "XE25I") | 
                     gp_clinical$read_2 %in% c("443e.", "4433"))
serum_LH$hormone <- "serum LH"

LH <- list(plasma_LH, serum_LH)

# Calculate age at measurement ----

cleaned <- lapply(LH, function (d) {
  
  d$event_dt <- as.Date(d$event_dt, "%d/%m/%Y")
  
  res <- subset(d, !is.na(d$event_dt))
  res <- subset(res, res$event_dt > res$dob)
  res <- subset(res, res$event_dt <= as.Date("2020-10-01"))
  
  res$age_years <- 
    age_calc(res$dob, res$event_dt, units = "years")
  
  return (res)
  
})

# Data metrics and cleaning ----

cleaned <- bind_rows(cleaned)

# Get break-up of LH measurements by sex and age
cleaned %>% group_by(hormone, sex) %>% count()

indivs <- cleaned %>% group_by(eid, sex, hormone) %>% count()
table(indivs[, c("sex", "hormone")])

plot_dat <- cleaned
plot_dat$value1 <- as.numeric(plot_dat$value1)
pdf("/well/lindgren/UKBIOBANK/samvida/hormones/LH/initial.pdf")
ggplot(plot_dat, aes(x = value1, color = sex)) +
  facet_wrap(~hormone) +
  geom_density(alpha = 0.25, na.rm = T) +
  labs(x = "Luteinising hormone", y = "Density")
dev.off()

# Add information on UKBIOBANK BMI ----

# Add BMI data from UKBIOBANK
ukbb_bmi <- pheno[, c("f.eid", 
                      "f.21001.0.0", "f.21001.1.0", "f.21001.2.0")]
ukbb_bmi$mean_UKBB_BMI <- rowMeans(ukbb_bmi[, 2:4], na.rm = T)

colnames(ukbb_bmi)[1] <- "eid"

cleaned <- merge(cleaned, ukbb_bmi[, c("eid", "mean_UKBB_BMI")],
                 by = "eid", all.x = T)

saveRDS(cleaned, "tmp.rds")

# Inconsistencies ----

# Convert value1/2/3 to numeric and only keep rows with numeric values

value_real <- apply(cleaned[, c("value1", "value2", "value3")], 1, function(x) {
  r <- as.numeric(x)
  r <- r[!is.na(as.numeric(x))]
  if (length(r) == 0) r <- NA
  return (r)
})

cleaned$value_real <- value_real
cleaned <- subset(cleaned, !is.na(cleaned$value_real))

colnames(cleaned)[16] <- "hormone_level"

# Remove duplicate rows
cleaned <- distinct(cleaned, 
                    eid, sex, mean_UKBB_BMI, data_provider, event_dt, 
                    age_years, hormone, hormone_level)

write.table(cleaned, "/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/LH_primary_care_annotated.txt",
            sep = "\t", quote = F, row.names = F)

## Cleaned stats ----

LH <- read.table("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/LH_primary_care_annotated.txt",
                 sep = "\t", header = T)

# Distribution of LH levels by sex
pdf("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/plots/sex_distribution.pdf")

ggplot(LH, aes(x = hormone_level, color = sex, fill = sex)) +
  facet_wrap(~hormone, nrow = 2) +
  geom_density(alpha = 0.25, na.rm = T) +
  labs(x = "LH Level", y = "Density") +
  xlim(c(0, 175))

ggplot(subset(LH, LH$sex == "F"), 
       aes(x = hormone, y = hormone_level)) +
  geom_jitter(size = 0.1, color = "grey") +
  geom_violin(fill = "#f8766d", trim = F) +
  ylim(c(0, 175)) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  labs(x = "Hormone", y = "Hormone Level", 
       title = "Females, Plasma LH n = 238, Serum LH n = 33,006")

ggplot(subset(LH, LH$sex == "M"), 
       aes(x = hormone, y = hormone_level)) +
  geom_jitter(size = 0.1, color = "grey") +
  geom_violin(fill = "#00bfc4", trim = F) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  labs(x = "Hormone", y = "Hormone Level", 
       title = "Males, Plasma LH n = 31, Serum LH n = 3,415")

dev.off()

# Distribution of LH values by age

breakpoints <- c(18, seq(30, 80, by = 10))
names <- c("(18-30]", "(30-40]", "(40-50]", "(50-60]", "(60-70]",
           "(70-80]")
LH$age_10_bin <- cut(LH$age_years, breaks = breakpoints, labels = names)

pdf("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/plots/age_distribution.pdf")

ggplot(subset(LH, LH$sex == "F"), 
       aes(x = age_10_bin, y = hormone_level)) +
  facet_wrap(~hormone, nrow = 2) +
  geom_violin(fill = "#f8766d") +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  ylim(c(0, 175)) +
  labs(x = "Age (years)", y = "Hormone Level", 
       title = "Females, Plasma LH n = 238, Serum LH n = 33,006")

ggplot(subset(LH, LH$sex == "M"), 
       aes(x = age_10_bin, y = hormone_level)) +
  facet_wrap(~hormone, nrow = 2, scales = "free_y") +
  geom_violin(fill = "#00bfc4", trim = F) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  labs(x = "Age (years)", y = "Hormone Level", 
       title = "Males, Plasma LH n = 31, Serum LH n = 3,415")

dev.off()

# Create categories for BMI (underweight, normal, overweight, etc.)
# based on WHO guidelines

breakpoints <- c(-Inf, 18.5, 25, 30, 35, 40, Inf)
names <- c("Underweight (< 18.5)", "Normal weight [18.5 - 25)", 
           "Pre-obesity [25 - 30)", "Obesity Class I [30 - 35)", 
           "Obesity Class II [35 - 40)", "Obesity Class III (>= 40)")
LH$BMI_class <- cut(LH$mean_UKBB_BMI, breaks = breakpoints, 
                    include.lowest = T, right = F,
                    labels = names)

pdf("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/plots/BMI_distribution.pdf")

ggplot(subset(LH, LH$sex == "F" & !is.na(LH$BMI_class)), 
       aes(x = BMI_class, y = hormone_level)) +
  facet_wrap(~hormone, nrow = 2) +
  geom_violin(fill = "#f8766d") +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  ylim(c(0, 175)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "UKBB BMI Class", y = "Hormone Level", 
       title = "Females, Plasma LH n = 238, Serum LH n = 33,006")

ggplot(subset(LH, LH$sex == "M" & !is.na(LH$BMI_class)), 
       aes(x = BMI_class, y = hormone_level)) +
  facet_wrap(~hormone, nrow = 2, scales = "free_y") +
  geom_violin(fill = "#00bfc4") +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "UKBB BMI Class", y = "Hormone Level", 
       title = "Males, Plasma LH n = 31, Serum LH n = 3,415")

dev.off()

######################################## STOP HERE FOR NOW #############

# Create categories for mean BMI (underweight, normal, overweight, etc.)
# based on WHO guidelines - how many individuals are in each category?

pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/BMI_classes.pdf")

bmi %>% distinct(eid, .keep_all = T) %>%
  group_by(sex, bmi_UKBIOBANK_class) %>% 
  count() %>%
  drop_na(bmi_UKBIOBANK_class) %>%
  ggplot(aes(x = bmi_UKBIOBANK_class, y = n, 
             color = sex, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "UKBIOBANK BMI Class", y = "Number of individuals")

bmi %>% distinct(eid, .keep_all = T) %>%
  group_by(sex, bmi_primarycare_class) %>% 
  count() %>%
  ggplot(aes(x = bmi_primarycare_class, y = n, 
             color = sex, fill = sex)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Primary Care BMI Class", y = "Number of individuals")

dev.off()

## Multiple measures ----

# Distribution of number of BMI events
indivs <- bmi %>% group_by(eid, sex) %>% count()
pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/number_events.pdf")

ggplot(indivs, aes(x = n, fill = sex, colour = sex)) +
  geom_histogram(alpha = 0.25, position = "identity") +
  labs(x = "Number of BMI measures", y = "Number of individuals")

ggplot(indivs, aes(x = n, fill = sex, colour = sex)) +
  geom_histogram(alpha = 0.25, position = "identity") +
  xlim(c(0, 50)) +
  labs(x = "Number of BMI measures", y = "Number of individuals")

dev.off()

# Number of measures by BMI class
pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/nmeasures_BMI_class.pdf")

tmp <- bmi %>% distinct(eid, .keep_all = T) 

ggplot(tmp, aes(x = bmi_primarycare_class, y = n, fill = sex)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Primary Care BMI Class", 
       y = "Number of BMI measures")

ggplot(tmp, aes(x = bmi_primarycare_class, y = n, fill = sex)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(c(0, 50)) +
  labs(x = "Primary Care BMI Class", 
       y = "Number of BMI measures")


dev.off()

# Count number of BMI measurements per individual
bmi <- bmi %>% group_by(eid) %>% add_tally()

# Create categories for number of measurements
breakpoints <- c(-Inf, 1, 5, 10, 15, Inf)
names <- c("1", "2-5", "6-10", "11-15", "16+")
bmi$nmeasures_bin <- cut(bmi$nn, breaks = breakpoints, labels = names)

# BMI values by number of measures
pdf("/well/lindgren/UKBIOBANK/samvida/BMI/plots/nmeasures_bin.pdf")

ggplot(bmi, aes(x = nmeasures_bin, y = primarycare_BMI, fill = sex)) +
  geom_boxplot() +
  ylim(c(10, 100)) +
  labs(x = "Number of BMI measures")

ggplot(bmi, aes(x = primarycare_BMI, fill = nmeasures_bin, 
                color = nmeasures_bin)) +
  geom_density(alpha = 0.25) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  xlim(c(10, 100))

dev.off()

# Is the relationship between BMI and number of measures due to age?

no_age_model <- lm(n ~ BMI, data = bmi)
summary(no_age_model) 
age_model <- lm(n ~ BMI + age_years, data = bmi)
summary(age_model) 

# Subset multiply measured individuals
bmi_mult <- subset(bmi, bmi$n > 1)

# Summary statistics for multiply measured individuals
bmi_mult <- bmi_mult %>% mutate(min_BMI = min(BMI), mean_BMI = mean(BMI), 
                                max_BMI = max(BMI), 
                                range_BMI = max_BMI - mean_BMI)
