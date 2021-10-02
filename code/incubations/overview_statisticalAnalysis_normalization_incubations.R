#### code/incubations/overview_statisticalAnalysis_normalization_incubations.R ####
# Benjamin D. Peterson

# Additional notes in the corresponding Obsidian document:
# Overview, statistical analysis, normalization of incubations

# Powerpoint with notes: results/overview_statisticalAnalysis_normalization_incubations

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")


#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                         sheet = "cleaned_incubation_data",) %>%
  mutate(coreID = fct_relevel(coreID, MG.order),
         matrixID = fct_relevel(matrixID, PW.order))


#### Plot MeHg production, faceted by sediment core (shows impact of porewater) ####
inc.Hg.data %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent) * 100,
            meth_spike_per_sd = sd(SMHG_201_percent),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  mutate(coreID = paste("Core source: ", coreID, sep = ""),
         coreID = fct_relevel(coreID, paste("Core source: ", MG.order, sep = ""))) %>%
  ggplot(aes(x = matrixID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = matrixID)) +
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  facet_wrap(~coreID, nrow = 2) +
  scale_fill_manual(values = color.vector) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black"))



#### Plot MeHg production, faceted by porewater source (shows impact of core) ####
inc.Hg.data %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent) * 100,
            meth_spike_per_sd = sd(SMHG_201_percent),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  mutate(matrixID = paste("Porewater source: ", matrixID, sep = ""),
         matrixID = fct_relevel(matrixID, paste("Porewater source: ", PW.order, sep = ""))) %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = coreID)) +
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  facet_wrap(~matrixID, nrow = 3) +
  scale_fill_manual(values = color.vector) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black"))


#### Test for normality ####
# Normality of overall data set
shapiro.test(inc.Hg.data$SMHG_201_percent)
hist(inc.Hg.data$SMHG_201_percent,
     breaks = 20)
# Not even close to normal.

# Log transformed data
shapiro.test(log(inc.Hg.data$SMHG_201_percent, 10))
hist(log(inc.Hg.data$SMHG_201_percent, 10),
     breaks = 20)
# Still not close.


#### Check validity of ANOVA as analytical tool ####

# Let's see if the residuals of an ANOVA are normally
# distributed, since that seems to be the main issue.

# Non-transformed data
twoWayAnova_inc <- aov(SMHG_201_percent ~ coreID * matrixID,
                       data = inc.Hg.data)
shapiro.test(twoWayAnova_inc$residuals)
hist(twoWayAnova_inc$residuals,
     breaks = 40)
qqnorm(twoWayAnova_inc$residuals)
# Looks bad...
rm(twoWayAnova_inc)


# Log-transformed data
twoWayAnova_inc_log <- aov(log(SMHG_201_percent, 10) ~ coreID * matrixID,
                           data = inc.Hg.data)
hist(twoWayAnova_inc_log$residuals,
     breaks = 20)
par(mfrow = c(1,2))
plot(density(twoWayAnova_inc_log$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
shapiro.test(twoWayAnova_inc_log$residuals)
# QQ-normal plot
qqnorm(twoWayAnova_inc_log$residuals)
qqline(twoWayAnova_inc_log$residuals)



#### Storming ahead with two-way ANOVA ####
summary(twoWayAnova_inc_log)



#### Plot MeHg production across cores on same axes, linked by porewater matrix ####
inc.Hg.data %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent) * 100,
            meth_spike_per_sd = sd(SMHG_201_percent),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             group = matrixID,
             col = matrixID)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = color.vector) +
  theme_bw() +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black"))



#### Plot MeHg production across porewater matrices on same axes, linked by core used ####
inc.Hg.data %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent) * 100,
            meth_spike_per_sd = sd(SMHG_201_percent),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  ggplot(aes(x = matrixID,
             y = meth_spike_per_mean,
             width = 0.8,
             group = coreID,
             col = coreID)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = color.vector) +
  theme_bw() +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black"))



#### Normalize data, both within groups of porewater and cores ####
inc.Hg.data.max.porewater <- inc.Hg.data %>%
  group_by(matrixID) %>%
  summarise(porewater_matrix_max_methylation = max(SMHG_201_percent))
inc.Hg.data.max.core <- inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(core_max_methylation = max(SMHG_201_percent))

inc.Hg.data.normalized <- full_join(inc.Hg.data,
                                    inc.Hg.data.max.porewater) %>%
  full_join(inc.Hg.data.max.core) %>%
  mutate(RMP_core = (SMHG_201_percent) / (porewater_matrix_max_methylation) * 100,
         RMP_porewater = (SMHG_201_percent) / (core_max_methylation) * 100)



#### Save out data ####
saveRDS(inc.Hg.data.normalized,
        "dataEdited/incubations/incubation_data_with_normalization.rds")
