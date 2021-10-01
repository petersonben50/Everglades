#### code/incubations ####
# Benjamin D. Peterson

# Obsidian notes here: Incubations - statistical analysis porewaters
# Results: results/RMP_porewater.pptx

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")



#### Read in plots ####
Hg.inc.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds")


#### Calculate core RMP mean and SE ####
coreID.data <- Hg.inc.data %>%
  group_by(coreID) %>%
  summarise(RMP_core_mean = mean(RMP_core),
            RMP_core_sd = sd(RMP_core),
            RMP_core_count = n(),
            RMP_core_se = RMP_core_sd / sqrt(RMP_core_count)) %>%
  ungroup() %>%
  rename(siteID = coreID) %>%
  select(siteID, RMP_core_mean, RMP_core_sd, RMP_core_se)


#### Calculate porewater RMP mean and SE ####
matrixID.data <- Hg.inc.data %>%
  group_by(matrixID) %>%
  summarise(RMP_porewater_mean = mean(RMP_porewater),
            RMP_porewater_sd = sd(RMP_porewater),
            RMP_porewater_count = n(),
            RMP_porewater_se = RMP_porewater_sd / sqrt(RMP_porewater_count)) %>%
  ungroup() %>%
  filter(matrixID %in% MG.order) %>%
  rename(siteID = matrixID) %>%
  select(siteID, RMP_porewater_mean, RMP_porewater_sd, RMP_porewater_se)


#### Calculate native MeHg production ####
native.MeHg.production.data <- Hg.inc.data %>%
  filter(as.character(coreID) == as.character(matrixID)) %>%
  group_by(coreID) %>%
  mutate(SMHG_201_percent = SMHG_201_percent * 100) %>%
  summarise(SMHG_201_mean = mean(SMHG_201_percent),
            SMHG_201_sd = sd(SMHG_201_percent),
            SMHG_201_count = n(),
            SMHG_201_se = SMHG_201_sd / sqrt(SMHG_201_count)) %>%
  ungroup() %>%
  rename(siteID = coreID) %>%
  select(siteID, SMHG_201_mean, SMHG_201_sd, SMHG_201_se)


#### Combine data ####
all.data <- coreID.data %>%
  full_join(matrixID.data) %>%
  full_join(native.MeHg.production.data)


#### Generate plot of MeHg production as function of two RMPs ####
all.data %>%
  ggplot(aes(x = RMP_core_mean,
             y = RMP_porewater_mean)) +
  geom_point(aes(size = SMHG_201_mean*2)) +
  scale_size_continuous(range = c(0.5, 6)) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlim(c(0, 100)) +
  ylab("Porewater RMP") +
  xlab("Sediment RMP") +
  labs(size = "MeHg production (%)") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.8, 0.8))

# Model response
all.data %>%
  ggplot(aes(x = RMP_core_mean,
             y = RMP_porewater_mean)) +
  geom_point(aes(size = (RMP_core_mean / 100)^2*(RMP_porewater_mean / 100)^2)) +
  scale_size_continuous(range = c(0.5, 6)) +
  theme_bw() +
  ylim(c(0, 100)) +
  xlim(c(0, 100)) +
  ylab("Porewater RMP") +
  xlab("Sediment RMP") +
  labs(size = "MeHg production (%)") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.8, 0.8))
# Holy shit this is so close to what the actual data looks like




#### 3D plot ####
scatterplot3d(data.frame(all.data$RMP_core_mean,
                         all.data$RMP_porewater_mean,
                         all.data$SMHG_201_mean),
              angle = 30,
              pch = 16,
              # size = 4,
              color = color.vector[all.data$siteID],
              type = "h",
              lab = c(4, 4),
              lab.z = 4,
              xlim = c(0, 100),
              ylim = c(0, 100),
              zlim = c(0, 3),
              xlab = "Core RMP",
              ylab = "Porewater RMP",
              zlab = "MeHg production under native conditions")