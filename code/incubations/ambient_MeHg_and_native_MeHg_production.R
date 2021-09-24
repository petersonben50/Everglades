#### code/incubations/ambient_MeHg_and_native_MeHg_production.R ####
# Benjamin D. Peterson


#### Get cleaned up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(patchwork)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                         sheet = "cleaned_incubation_data",) %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) 


#### Read in porewater data, order it, and generate percent MeHg ####
porewater.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                            sheet = "porewater_data") %>%
  filter(siteID %in% MG.order) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  mutate(MeHg_porewater = as.numeric(FMHG),
         HgT_porewater = as.numeric(FTHG),
         percent_MeHg_porewater = MeHg_porewater / HgT_porewater * 100) %>%
  select(siteID, MeHg_porewater, HgT_porewater, percent_MeHg_porewater)



#### Calculate ambient MeHg in sediment from incubation cores ####
amb.Hg.sediment <- inc.Hg.data %>%
  group_by(coreID) %>%
  mutate(SMHG_201_percent = SMHG_201_percent * 100) %>%
  summarise(SMHG_amb_mean = mean(SMHG_amb),
            SMHG_amb_sd = sd(SMHG_amb),
            SMHG_amb_count = n(),
            SMHG_amb_se = SMHG_amb_sd / sqrt(SMHG_amb_count),
            STHG_amb_mean = mean(STHG_amb),
            STHG_amb_sd = sd(STHG_amb),
            STHG_amb_count = n(),
            STHG_amb_se = STHG_amb_sd / sqrt(STHG_amb_count),
            SMHG_amb_percent_mean = mean(SMHG_amb_percent),
            SMHG_amb_percent_sd = sd(SMHG_amb_percent),
            SMHG_amb_percent_count = n(),
            SMHG_amb_percent_se = SMHG_amb_percent_sd / sqrt(SMHG_amb_percent_count)) %>%
  select(coreID, SMHG_amb_mean, SMHG_amb_se,
         STHG_amb_mean, STHG_amb_se,
         SMHG_amb_percent_mean, SMHG_amb_percent_se)
amb.Hg.sediment



#### Generate plot of ambient MeHg in soil ####
MeHg.soil.plot <- amb.Hg.sediment %>%
  ggplot(aes(x = coreID,
             y = SMHG_amb_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = SMHG_amb_mean - SMHG_amb_se,
                    ymax = SMHG_amb_mean + SMHG_amb_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(title = "MeHg concentrations in sediment",
       y = "MeHg (ng/g)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Generate plot of percent ambient MeHg ####
MeHgPercent.soil.plot <- amb.Hg.sediment %>%
  ggplot(aes(x = coreID,
             y = SMHG_amb_percent_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = SMHG_amb_percent_mean - SMHG_amb_percent_se,
                    ymax = SMHG_amb_percent_mean + SMHG_amb_percent_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(title = "Percent MeHg in sediment",
       y = "Percent MeHg") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Generate plot of HgT levels ####
HgT.soil.plot <- amb.Hg.sediment %>%
  ggplot(aes(x = coreID,
             y = STHG_amb_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = STHG_amb_mean - STHG_amb_se,
                    ymax = STHG_amb_mean + STHG_amb_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(title = "HgT in sediment",
       y = "HgT (ng/g)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Stack plots together ####
MeHg.soil.plot / HgT.soil.plot / MeHgPercent.soil.plot
rm(MeHg.soil.plot, MeHgPercent.soil.plot, HgT.soil.plot)




#### Porewater MeHg concentration plot ####
MeHg.porewater.plot <- porewater.data %>%
  ggplot(aes(x = siteID,
             y = MeHg_porewater,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector) +
  labs(title = "MeHg in porewater",
       y = "MeHg (ng/L)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Porewater HgT concentration plot ####
HgT.porewater.plot <- porewater.data %>%
  ggplot(aes(x = siteID,
             y = HgT_porewater,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector) +
  labs(title = "HgT in porewater",
       y = "HgT (ng/L)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Porewater MeHg percentage plot ####
MeHg.percent.porewater.plot <- porewater.data %>%
  ggplot(aes(x = siteID,
             y = percent_MeHg_porewater,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector) +
  labs(title = "Percent MeHg in porewater",
       y = "MeHg (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Porewater plots stacked ####
MeHg.porewater.plot / HgT.porewater.plot / MeHg.percent.porewater.plot
rm(MeHg.porewater.plot, HgT.porewater.plot, MeHg.percent.porewater.plot)



#### Production MeHg using native porewaters calculations ####
native.porewater.production <- inc.Hg.data %>%
  filter(coreID == matrixID) %>%
  group_by(coreID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent * 100),
            meth_spike_per_sd = sd(SMHG_201_percent * 100),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count)) %>%
  select(coreID, meth_spike_per_mean, meth_spike_per_se)


#### Plot MeHg production using native cores ####
native.porewater.production %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(y = "Methylated spike (%)",
       title = "Hg methylation using ambient porewater") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))



#### Native MeHg production vs. MeHg in porewater ####
# Generate data frame for MeHg percent in porewater against MeHg production
MeHg.prod.and.porewater <- native.porewater.production %>%
  left_join(porewater.data %>% rename(coreID = siteID))
# Make plot
MeHg.prod.vs.porewater.MeHg.levels <- MeHg.prod.and.porewater %>%
  ggplot(aes(x = MeHg_porewater,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector) +
  labs(x = "MeHg (ng/L)",
       y = "Methylated spike (%)",
       title = "Native Hg methylation vs. porewater MeHg levels") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
# Make plot
MeHg.prod.vs.porewater.MeHg.percent <- MeHg.prod.and.porewater %>%
  ggplot(aes(x = percent_MeHg_porewater,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector) +
  labs(x = "MeHg (%)",
       y = "Methylated spike (%)",
       title = "Native Hg methylation vs. porewater MeHg percent") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))



#### Native MeHg production vs. MeHg in sediment ####
# MeHg percent in sediment against MeHg production
MeHg.prod.and.sediment <- native.porewater.production %>%
  left_join(amb.Hg.sediment)


# Plot of MeHg production vs. MeHg levels in seds
MeHg.prod.vs.sediment.MeHg.levels <- MeHg.prod.and.sediment %>% 
  ggplot(aes(x = SMHG_amb_mean,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector) +
  labs(y = "Methylated spike (%)",
       x = "MeHg levels (ng/g)",
       title = "Native Hg methylation vs. sediment MeHg levels") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
# Plot of MeHg production vs. MeHg percent in seds
MeHg.prod.vs.sediment.MeHg.percent <- MeHg.prod.and.sediment %>% 
  ggplot(aes(x = SMHG_amb_percent_mean,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector) +
  labs(y = "Methylated spike (%)",
       x = "MeHg (%)",
       title = "Native Hg methylation vs. sediment MeHg percent") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))


#### Plot figures together ####
(MeHg.prod.vs.porewater.MeHg.levels + MeHg.prod.vs.porewater.MeHg.percent) / (MeHg.prod.vs.sediment.MeHg.levels + MeHg.prod.vs.sediment.MeHg.percent)


#### Run linear regressions of all combos ####

# Native MeHg production vs sediment MeHg concentration
summary(lm(meth_spike_per_mean ~ MeHg_porewater,
           data = MeHg.prod.and.porewater))
# No linear regression

# Native MeHg production vs sediment MeHg percent
summary(lm(meth_spike_per_mean ~ percent_MeHg_porewater,
           data = MeHg.prod.and.porewater))
# No linear regression

# Native MeHg production vs sediment MeHg concentration
summary(lm(meth_spike_per_mean ~ SMHG_amb_mean,
           data = MeHg.prod.and.sediment))
# No linear regression

# Native MeHg production vs sediment MeHg percent
summary(lm(meth_spike_per_mean ~ SMHG_amb_percent_mean,
           data = MeHg.prod.and.sediment))
# No linear regression
