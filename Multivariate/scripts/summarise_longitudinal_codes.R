# Author: Samvida S. Venkatesh
# Date: 15/01/2021

library(tidyverse)
theme_set(theme_bw())

# Generate merged codes ----

v2v3 <- read.table("/well/lindgren/UKBIOBANK/samvida/v2v3_merge.txt",
                   sep = "\t", header = T, stringsAsFactors = F,
                   quote = "", fill = F, comment.char = "~")

v2v3_merge <- v2v3[, c("READV2_CODE", "READV3_CODE")]
v2v3_merge$DESC <- ifelse(v2v3$TERMV2_DESC == v2v3$TERMV3_DESC, 
                          v2v3$TERMV2_DESC,
                          paste0(v2v3$TERMV2_DESC, ";", v2v3$TERMV3_DESC))
v2v3_merge$unique_code <- seq(1:dim(v2v3_merge)[1])

v2v3_merge <- v2v3_merge %>% distinct(READV2_CODE, READV3_CODE, .keep_all = T)

write.table(v2v3_merge, "C:/Users/samvida/Documents/Lindgren Group/Resources/UKBIOBANK/Primary Care/merged_v2v3_codes.txt",
            sep = "\t", quote = F, row.names = F)

# Grab longitudinal multivariate with merged codes ----

v2v3_merge <- read.table("/well/lindgren/UKBIOBANK/samvida/merged_v2v3_codes.txt",
                         sep = "\t", header = T, stringsAsFactors = F,
                         quote = "", fill = F, comment.char = "~")

gp_clinical <- readRDS("/well/lindgren/UKBIOBANK/samvida/gp_clinical_annotated_sex_dob.rds")

# Add whether there is a quantitative measurement for the trait
gp_clinical$quant <- apply(gp_clinical[, c("value1", "value2", "value3")], 1, 
                           function (x) { any(grepl("[0-9].", x)) } )

# Add unique code from read 2
gp_clin_new <- gp_clinical[, c("eid", "event_dt", "read_2", 
                               "read_3", "quant")]
match_ind <- match(gp_clin_new$read_2, v2v3_merge$READV2_CODE)
gp_clin_new$unique_code_v2 <- ifelse(is.na(match_ind), gp_clin_new$read_2, 
                                     v2v3_merge$unique_code[match_ind])

# Add unique code from read 3
match_ind <- match(gp_clin_new$read_3, v2v3_merge$READV3_CODE)
gp_clin_new$unique_code_v3 <- ifelse(is.na(match_ind), gp_clin_new$read_3, 
                                     v2v3_merge$unique_code[match_ind])
# Merge unique codes
gp_clin_new$unique_code <- 
  ifelse(gp_clin_new$unique_code_v2 == gp_clin_new$unique_code_v3,
         gp_clin_new$unique_code_v2,
         paste0(gp_clin_new$unique_code_v2, gp_clin_new$unique_code_v3))

summary <- gp_clin_new %>% count(unique_code, quant, eid) 

# Only keep codes for which we have data at >1 timepoint on >=1 individual
keep <- summary %>% group_by(unique_code, quant) %>% 
  summarise(n_indivs = n(), max_measures = max(n))
keep <- subset(keep, keep$n_indivs >= 1 & keep$max_measures > 1)$unique_code

summary$code_quant <- ifelse(summary$quant, paste0(summary$unique_code, "_q"),
                             summary$unique_code)

summary <- subset(summary, summary$unique_code %in% keep) %>% 
  split(.$code_quant)
summary <- lapply(summary, function (v) {
  return (setNames(v$n, v$eid))
})

saveRDS(summary, "/well/lindgren/UKBIOBANK/samvida/multivariate/unique_codes_summary.rds")

# Examine with summary statistics ----

summary <- readRDS("/well/lindgren/UKBIOBANK/samvida/multivariate/unique_codes_summary.rds")
unique_code_lkp <- read.table("/well/lindgren/UKBIOBANK/samvida/merged_v2v3_codes.txt",
                              sep = "\t", header = T, stringsAsFactors = F,
                              quote = "", fill = F, comment.char = "~")

dat <- data.frame(unique_code = unlist(names(summary)))
dat <- separate(dat, unique_code, into = c("unique_code", "quant"), sep = "_")
dat$quant <- !is.na(dat$quant)

# Calculate 
dat$n_indivs <- unlist(lapply(summary, function (x) length(x)))
dat$mean_measures <- unlist(lapply(summary, function (x) mean(x)))
dat$median_measures <- unlist(lapply(summary, function (x) median(x)))
dat$max_measures <- unlist(lapply(summary, function (x) max(x)))

N_codes <- dim(dat)[1]
# add description column
desc_ind <- match(dat$unique_code, unique_code_lkp$unique_code)
no_match <- which(is.na(desc_ind))
# check if it matches a V2 code
v2_match <- match(dat$unique_code[no_match],
                  unique_code_lkp$READV2_CODE)
desc_ind[no_match] <- v2_match
no_match <- which(is.na(desc_ind))
# check if it matches a v3 code
desc_ind[no_match] <- match(dat$unique_code[no_match],
                            unique_code_lkp$READV3_CODE)

dat$term_description <- ifelse(is.na(desc_ind), NA, 
                               unique_code_lkp$DESC[desc_ind])

write.table(dat, "summary_stats.txt", sep = "\t", quote = F,
            row.names = F)

## Scatter plots of n_indivs vs mean, median, max obs per trait ----

# Long format for plots
dat <- pivot_longer(dat, cols = ends_with("measures"),
                    names_to = "stat",
                    values_to = "value")
dat$label_desc <- ifelse(dat$quant, paste0(dat$term_description, ", quant"),
                         dat$term_description)

# Label interesting points in mean and median measures 
sub_dat <- subset(dat, dat$stat != "max_measures")
sub_dat$label <- sub_dat$n_indivs >= 10000 & sub_dat$value >= 10

pdf("scatter_plots.pdf", onefile = T)

ggplot(subset(sub_dat, sub_dat$stat == "mean_measures"), 
       aes(x = n_indivs, y = value)) +
  facet_wrap(~stat) +
  geom_point() +
  geom_text_repel(data = . %>% filter(label), 
            aes(label = label_desc)) + 
  labs(x = "# individuals", 
       y = "Mean # observations per individual",
       title = paste0("Combined read codes, N = ", N_codes))

dev.off()
