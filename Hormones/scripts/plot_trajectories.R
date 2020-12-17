library(tidyverse)
library(lubridate)
theme_set(theme_bw())
library(RColorBrewer)

LH <- read.table("LH_primary_care_annotated.txt", 
                 sep = "\t", header = T, stringsAsFactors = F)

# Correlation between serum and plasma LH ----

serum_ids <- unique(LH$eid[LH$hormone == "serum LH"])
plasma_ids <- unique(LH$eid[LH$hormone == "plasma LH"])
both_ids <- intersect(serum_ids, plasma_ids)

both <- LH[LH$eid %in% both_ids, ]
both <- pivot_wider(both, id_cols = c("eid", "sex", "event_dt"), 
                    names_from = "hormone", values_from = "hormone_level")

both <- both[which(both$`plasma LH` != "NULL" & 
                     both$`serum LH` != "NULL"), ]

# Only 20 individuals with both, all of them have 0 in one or the other!

# Longitudinal summaries ----

indivs <- LH %>% group_by(eid, sex, hormone) %>% 
  summarise(n_obs = n())

# Group individuals by number of LH measures (1, 2-5, 6-10, 11+)
# Create categories for number of measurements
breakpoints <- c(-Inf, 1, 5, 10, Inf)
names <- c("1", "2-5", "6-10", "11+")
indivs$nmeasures_bin <- cut(indivs$n_obs, breaks = breakpoints, labels = names)

pdf("plots/longitudinal_summaries.pdf")

# Histogram of number of LH measures
ggplot(indivs, aes(x = n_obs, fill = sex, colour = sex)) +
  facet_wrap(~hormone, ncol = 2, scales = "free") +
  geom_histogram(alpha = 0.25, position = "identity") +
  labs(x = "Number of LH measures", y = "Number of individuals")

dev.off()

# Randomly sample trajectories in each nobs_bin ----

# Number of LH measurements
LH <- LH %>% group_by(eid, hormone) %>% mutate(n_obs = n())

# Group by number of measurements
breakpoints <- c(-Inf, 1, 5, 10, Inf)
names <- c("1", "2-5", "6-10", "11+")
LH$nobs_bin <- cut(LH$n_obs, breaks = breakpoints, labels = names)

LHlong <- subset(LH, LH$nobs_bin != "1")

pids <- LHlong %>% distinct(eid, sex, hormone, nobs_bin) %>%
  group_by(sex, hormone, nobs_bin) %>%
  sample_n(size = min(n(), 10))

pdata <- LHlong[LHlong$eid %in% pids$eid, ]

pdf("plots/trajectories.pdf", onefile = T)

# Females
ggplot(subset(pdata, pdata$sex == "F"), 
       aes(x = age_years, y = hormone_level, group = eid)) +
  facet_wrap(~nobs_bin+hormone, nrow = 2, scales = "free_y") +
  geom_point(col = "#F8766D") +
  geom_line(col = "#F8766D", size = 0.7) +
  labs(x = "Age (years)", y = "LH Level")

# Males
ggplot(subset(pdata, pdata$sex == "M"), 
       aes(x = age_years, y = hormone_level, group = eid)) +
  facet_wrap(~nobs_bin+hormone, nrow = 2, scales = "free_y") +
  geom_point(col = "#00BFC4") +
  geom_line(col = "#00BFC4", size = 0.7) +
  labs(x = "Age (years)", y = "LH Level")

dev.off()

# Trajectories plotted from first observation ----

LHlong <- subset(LH, LH$nobs_bin != "1")

LHlong <- LHlong %>% group_by(eid, hormone) %>%
  arrange(age_years, .by_group = T) %>%
  mutate(interval_weeks = interval(first(event_dt), event_dt) %/% weeks(1))

pids <- LHlong %>% distinct(eid, sex, hormone, nobs_bin) %>%
  group_by(sex, hormone, nobs_bin) %>%
  sample_n(size = min(n(), 10))

pdata <- LHlong[LHlong$eid %in% pids$eid, ]

pdf("plots/trajectories_weeks_from_initial.pdf", onefile = T)

# Females
ggplot(subset(pdata, pdata$sex == "F"), 
       aes(x = interval_weeks, y = hormone_level, group = eid)) +
  facet_wrap(~nobs_bin+hormone, nrow = 2, scales = "free") +
  geom_point(col = "#F8766D") +
  geom_line(col = "#F8766D", size = 0.7) +
  labs(x = "Weeks from first measurement", y = "LH Level")

# Males
ggplot(subset(pdata, pdata$sex == "M"), 
       aes(x = interval_weeks, y = hormone_level, group = eid)) +
  facet_wrap(~nobs_bin+hormone, nrow = 2, scales = "free") +
  geom_point(col = "#00BFC4") +
  geom_line(col = "#00BFC4", size = 0.7) +
  labs(x = "Weeks from first measurement", y = "LH Level")

dev.off()

