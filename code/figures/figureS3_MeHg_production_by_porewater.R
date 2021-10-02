#### code/figures/figureS3_MeHg_production_by_porewater.R ####
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


#### Read in incubation data ####
inc.Hg.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds")


#### Order by average RMP ####
inc.Hg.data.order <- inc.Hg.data %>%
  group_by(matrixID) %>%
  summarise(RMP_porewater = mean(RMP_porewater)) %>%
  arrange(desc(RMP_porewater)) %>%
  select(matrixID) %>%
  mutate(matrixID = as.character(matrixID)) %>%
  unlist()
inc.Hg.data <- inc.Hg.data %>%
  mutate(matrixID = fct_relevel(matrixID,
                                inc.Hg.data.order))
color.vector <- color.vector[inc.Hg.data.order]

#### Plot MeHg production, faceted by sediment core, ordered by RMP ####
inc.data.plot.by.porewater <- inc.Hg.data %>%
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
  theme(axis.text.y = element_text(colour="black",
                                   size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())


#### Save out plot ####
pdf("results/figures/S3_MeHg_production_ordered_by_RMP.pdf",
    width = 12,
    height = 6)
inc.data.plot.by.porewater
dev.off()
