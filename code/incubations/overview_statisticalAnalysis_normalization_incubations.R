#### code/incubations/overview_statisticalAnalysis_normalization_incubations.R ####
# Benjamin D. Peterson

# Additional notes in the corresponding Obsidian document:
# Overview, statistical analysis, normalization of incubations

# Powerpoint with notes: results/overview_statisticalAnalysis_normalization_incubations

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(AICcmodavg)
library(ggpubr)
library(lme4)
library(readxl)
library(tidyverse)
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




#### Check validity of ANOVA as analytical tool ####

# Let's see if the residuals of an ANOVA are normally
# distributed, since that seems to be the main issue.

#### Non-transformed data ####
twoWayAnova_inc <- aov(SMHG_201_percent ~ coreID * matrixID,
                       data = inc.Hg.data)
shapiro.test(twoWayAnova_inc$residuals)
hist(twoWayAnova_inc$residuals,
     breaks = 40)
qqnorm(twoWayAnova_inc$residuals)
# Looks bad...
rm(twoWayAnova_inc)


#### Log-transformed data ####
twoWayAnova_inc_log <- aov(log(SMHG_201_percent, 10) ~ coreID * matrixID,
                           data = inc.Hg.data)
hist(twoWayAnova_inc_log$residuals,
     breaks = 20)
par(mfrow = c(1,2))
plot(density(twoWayAnova_inc_log$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot, test for normality of residuals
qqnorm(twoWayAnova_inc_log$residuals)
qqline(twoWayAnova_inc_log$residuals)

par(mfrow=c(2,2))
plot(twoWayAnova_inc_log)
# Heavy tailed.


#### Folded-root transformation ####
twoWayAnova_inc_fold <- aov((sqrt(SMHG_201_percent) - sqrt(100-SMHG_201_percent)) ~ coreID * matrixID,
                           data = inc.Hg.data)
hist(twoWayAnova_inc_fold$residuals,
     breaks = 20)
par(mfrow = c(2,2))
plot(density(twoWayAnova_inc_fold$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot, test for normality of residuals
qqnorm(twoWayAnova_inc_fold$residuals)
qqline(twoWayAnova_inc_fold$residuals)
par(mfrow=c(2,2))
plot(twoWayAnova_inc_fold)
# Still pretty heavy tailed.


#### Convert to rankings ####
inc.Hg.data.ranked <- inc.Hg.data %>%
  arrange(desc(SMHG_201_percent)) %>%
  mutate(rank = 1:length(SMHG_201_percent))

ggarrange(inc.Hg.data.ranked %>%
            ggplot(aes(y = rank,
                       x = coreID)) +
            geom_point() + theme_bw(),
          inc.Hg.data.ranked %>%
            ggplot(aes(y = rank,
                       x = matrixID)) +
            geom_point() + theme_bw())

twoWayAnova_inc_rank <- aov(rank ~ coreID * matrixID,
                            data = inc.Hg.data.ranked)
hist(twoWayAnova_inc_rank$residuals,
     breaks = 20)
par(mfrow = c(2,2))
plot(density(twoWayAnova_inc_rank$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot, test for normality of residuals
qqnorm(twoWayAnova_inc_rank$residuals)
qqline(twoWayAnova_inc_rank$residuals)

# Still pretty heavy tailed.
# Okay, no ANOVA. Skip to permutation test.


#### ANOVA with permutation tests ####
source("code/incubations/USP.r")
source("code/incubations/synchro_summary.R")
x <- inc.Hg.data %>%
  select(coreID, matrixID) %>%
  mutate(coreID = as.character(coreID),
         matrixID = as.character(matrixID)) %>%
  as.data.frame()
y <- inc.Hg.data %>%
  mutate(SMHG_201_percent_log = log(SMHG_201_percent, 10)) %>%
  select(SMHG_201_percent_log) %>%
  unlist()
#t <- USP(y, x, C = 1000000)
# saveRDS(object = t,
#         file = "dataEdited/incubations/synchronized_permutation_test_results.rds")
# Took a long time to run 1 million permutations, so saved out the R object
# for easy viewing later.
t <- readRDS("dataEdited/incubations/synchronized_permutation_test_results.rds")
synchro.summary(t)



#### Generate interaction plot ####
means.of.data <- inc.Hg.data %>%
  mutate(SMHG_201_percent_log = log(SMHG_201_percent, 10)) %>%
  group_by(coreID, matrixID) %>%
  summarise(SMHG_201_percent_log_mean = mean(SMHG_201_percent_log)) %>%
  select(coreID, matrixID, SMHG_201_percent_log_mean)
  
interaction.plot <- inc.Hg.data %>%
  mutate(SMHG_201_percent_log = log(SMHG_201_percent, 10)) %>%
  left_join(means.of.data) %>%
  ggplot(aes(x = coreID,
             width = 0.8,
             group = matrixID,
             col = matrixID)) +
  geom_point(aes(y = SMHG_201_percent_log)) +
  geom_line(aes(y = SMHG_201_percent_log_mean)) +
  scale_color_manual(values = color.vector) +
  theme_bw() +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black"))



#### Use AIC to test for interaction ####
Cand.mod <- list()
Cand.mod[["without_interaction"]] <- lm(log(SMHG_201_percent, 10) ~ coreID + matrixID,
                                        data = inc.Hg.data)
Cand.mod[["with_interaction"]] <- lm(log(SMHG_201_percent, 10) ~ coreID * matrixID,
                                     data = inc.Hg.data)
Cand.mod[["coreID"]] <- lm(log(SMHG_201_percent, 10) ~ coreID,
                                     data = inc.Hg.data)
Cand.mod[["matrixID"]] <- lm(log(SMHG_201_percent, 10) ~ matrixID,
                           data = inc.Hg.data)
aictab(cand.set = Cand.mod,modnames = names(Cand.mod))


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
