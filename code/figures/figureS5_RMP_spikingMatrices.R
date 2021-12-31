#### code/figures/figureS4_RMP_spikingMatrices.R ####
# Benjamin D. Peterson


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


#### Generate plot of RMP ####
RMP.boxplot <- inc.Hg.data %>%
  ggplot(aes(x = matrixID,
             y = RMP_porewater)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = coreID,
                  shape = coreID),
              width = 0.08) +
  scale_shape_manual(values = point.vector,
                     name = "Pore water matrix") +
  scale_color_manual(values = color.vector,
                     name = "Pore water matrix") +
  labs(y = "Pore water RMP (%)",
       x = "Pore water ID") +
  theme_bw() +
  theme(axis.text.y = element_text(colour="black",
                                   size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(colour="black",
                                   size = 11),
        axis.title.x = element_text(size = 14),
        legend.position = "none")
RMP.boxplot

#### Save out image ####
pdf("results/figures/S4_RMP_poreWater.pdf",
    width = 10,
    height = 6)
RMP.boxplot
dev.off()
