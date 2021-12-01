# Author: Samvida S. Venkatesh
# Date: 21/08/20

library(tidyr)
library(dplyr)
library(eeptools)
library(lubridate)
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

### Distribution of LH levels by sex ----

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

### Distribution of LH values by age ----

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

### Distribution by BMI ----

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

# Number of hormone measurements
LH <- LH %>% group_by(eid, hormone) %>% mutate(n_obs = n())

breakpoints <- c(-Inf, 1, 3, 9, Inf)
names <- c("1", "2-3", "4-9", "10+")
LH$nmeasures_bin <- cut(LH$n_obs, breaks = breakpoints, labels = names)

pdf("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/plots/nmeasures_distribution.pdf")

# LH distribution faceted by number of LH measures
ggplot(subset(LH, LH$sex == "F"), 
       aes(x = nmeasures_bin, y = hormone_level)) +
  facet_wrap(~hormone, nrow = 2) +
  geom_violin(fill = "#f8766d") +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  ylim(c(0, 175)) +
  labs(x = "Number of LH measures", y = "Hormone Level", 
       title = "Females, Plasma LH n = 238, Serum LH n = 33,006")

ggplot(subset(LH, LH$sex == "M"), 
       aes(x = nmeasures_bin, y = hormone_level)) +
  facet_wrap(~hormone, nrow = 2, scales = "free_y") +
  geom_violin(fill = "#00bfc4") +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.1,
               fill = "white") +
  labs(x = "Number of LH measures", y = "Hormone Level", 
       title = "Males, Plasma LH n = 31, Serum LH n = 3,415")

# BMI distribution faceted by number of LH measures
indivs <- distinct(LH, eid, sex, mean_UKBB_BMI, nmeasures_bin) 

ggplot(indivs, 
       aes(x = mean_UKBB_BMI, col = nmeasures_bin, fill = nmeasures_bin)) +
  facet_wrap(~sex, nrow = 2) +
  geom_density(alpha = 0.25, position = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "UKBIOBANK BMI", y = "Density")

dev.off()

### Number of tests in each age bin and BMI class etc. ----

pdf("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/plots/ntests.pdf")

# Age
ntests <- LH %>% group_by(sex, hormone, age_10_bin) %>% count()
ggplot(ntests, aes(x = age_10_bin, y = n, fill = sex, col = sex)) +
  facet_wrap(~hormone+sex, nrow = 2, scales = "free_y") +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Age (years)", y = "# tests performed")

# BMI
ntests <- LH %>% group_by(sex, hormone, BMI_class) %>% count()
ggplot(ntests, aes(x = BMI_class, y = n, fill = sex, col = sex)) +
  facet_wrap(~hormone+sex, nrow = 2, scales = "free_y") +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "BMI class", y = "# tests performed")

dev.off()

# Remove noise ----

# Calculate outliers
cutoffs <- LH %>% group_by(sex, hormone) %>% summarise(n = n(), 
                                                       five_indivs_cutoff = 5/n)
thresholds <- LH %>% group_by(sex, hormone) %>% 
  summarise(q1 = quantile(hormone_level, 0.9)[[1]], 
            q2 = quantile(hormone_level, 0.99)[[1]],
            q3 = quantile(hormone_level, 0.999)[[1]],
            q4 = quantile(hormone_level, 0.9999)[[1]])

# Manual inspection
pF <- subset(LH, LH$hormone == "plasma LH" & LH$sex == "F")
pF <- pF[order(pF$hormone_level, decreasing = T), ]
# Values over 150 are obvious errors - discard
pM <- subset(LH, LH$hormone == "plasma LH" & LH$sex == "M")
pM <- pM[order(pM$hormone_level, decreasing = T), ]
# No noise

sF <- subset(LH, LH$hormone == "serum LH" & LH$sex == "F")
sF <- sF[order(sF$hormone_level, decreasing = T), ]
# Values over 150 are obvious errors - discard
sM <- subset(LH, LH$hormone == "serum LH" & LH$sex == "M")
sM <- sM[order(sM$hormone_level, decreasing = T), ]
# No noise

cleaned <- subset(LH, LH$hormone_level < 150)
write.table(cleaned, 
            "/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH/LH_primary_care_annotated.txt", 
            sep = "\t", quote = F, row.names = F)