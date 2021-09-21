#### code/incubations/2019_inc_amb.R ####
# Benjamin D. Peterson

#### Get cleaned up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(patchwork)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Generate vector of correct order of samples along sulfate gradient ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")


#### Prep color vector ####
color.vector.core <- cb.translator[1:6]
names(color.vector.core) <- MG.order


#### Read in incubation data and order it ####
# inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis of core experiments_v2.xlsx",
#                          sheet = "cleaned_incubation_data") %>%
#   mutate(coreID = fct_relevel(coreID, MG.order))

inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                         sheet = "cleaned_incubation_data",) %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) 


#### Generate plot of ambient MeHg ####
# pdf("results/2019_incubations/ambient_MeHg.pdf",
#     width = 6,
#     height = 3)
inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(MeHg_amb_mean = mean(SMHG_amb),
            MeHg_amb_sd = sd(SMHG_amb),
            MeHg_amb_count = n(),
            MeHg_amb_se = MeHg_amb_sd / sqrt(MeHg_amb_count)) %>%
  ggplot(aes(x = coreID,
             y = MeHg_amb_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = MeHg_amb_mean - MeHg_amb_se,
                    ymax = MeHg_amb_mean + MeHg_amb_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector.core) +
  labs(title = "Ambient MeHg in cores",
       y = "Percent MeHg") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
# dev.off()



#### Generate plot of percent ambient MeHg ####
# pdf("results/2019_incubations/ambient_MeHg_fraction.pdf",
#     width = 6,
#     height = 3)
inc.Hg.data %>%
  mutate(SMHG_201_percent = SMHG_201_percent * 100) %>%
  group_by(coreID) %>%
  summarise(SMHG_201_percent_mean = mean(SMHG_201_percent),
            SMHG_201_percent_sd = sd(SMHG_201_percent),
            SMHG_201_percent_count = n(),
            SMHG_201_percent_se = SMHG_201_percent_sd / sqrt(SMHG_201_percent_count)) %>%
  ggplot(aes(x = coreID,
             y = SMHG_201_percent_mean,
             width = 0.8,
             fill = coreID)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = SMHG_201_percent_mean - SMHG_201_percent_se,
                    ymax = SMHG_201_percent_mean + SMHG_201_percent_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector.core) +
  labs(title = "Percent ambient MeHg in cores",
       y = "Percent MeHg") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
# dev.off()



#### Save out ambient MeHg data ####
# amb.Hg.data %>%
#   saveRDS("dataEdited/2019_incubations/ambient_MeHg.rds")


#### Plot out MeHg production using native porewaters ####

inc.Hg.data %>%
  filter(coreID == matrixID) %>%
  group_by(coreID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent * 100),
            meth_spike_per_sd = sd(SMHG_201_percent * 100),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count)) %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector.core) +
  labs(y = "Methylated spike (%)",
       title = "Hg methylation using ambient porewater") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))

#### Save out methylation data under native conditions ####
# inc.Hg.data %>%
#   saveRDS("dataEdited/2019_incubations/native_methylation.rds")
# 




#### Read in porewater data ####
PW.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                     sheet = "porewater_data",) %>%
  filter(siteID %in% MG.order) %>%
  mutate(siteID = fct_relevel(siteID, MG.order),
         FTHG = as.numeric(FTHG),
         FMHG = as.numeric(FMHG)) %>%
  mutate(FMHG_per = FMHG / FTHG)



#### MeHg in porewater ####

# MeHg percent in porewater
PW.MeHg <- PW.data %>%
  ggplot(aes(x = siteID,
             y = FMHG,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector.core) +
  labs(title = "Ambient MeHg levels in cores",
       y = "MeHg (ng/L)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))

# MeHg percent in porewater
PW.MeHg.per <- PW.data %>%
  ggplot(aes(x = siteID,
             y = FMHG_per,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector.core) +
  labs(title = "Ambient MeHg percent in cores",
       y = "MeHg (ng/L)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))

# Plot them
PW.MeHg / PW.MeHg.per



#### Plot MeHg production with native porewater against MeHg in porewater ####

# MeHg percent in porewater
inc.Hg.data %>%
  filter(coreID == matrixID) %>%
  group_by(coreID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent * 100),
            meth_spike_per_sd = sd(SMHG_201_percent * 100),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count)) %>%
  left_join(PW.data %>% rename(coreID = siteID)) %>%
  ggplot(aes(x = FMHG,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector.core) +
  labs(y = "Methylated spike (%)",
       title = "Hg methylation using ambient porewater") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
