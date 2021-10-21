#### code/figures/figureS6_alpha_diversity.R ####
# Benjamin D. Peterson


#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
# source("code/setup_PW_core_order_color_points.R")



richness.plot <- readRDS("results/metagenomes/assembly/hgcA/richness.rds")
E.plot <- readRDS("results/metagenomes/assembly/hgcA/evenness.rds")



pdf("results/figures/S6_alpha_diversity.pdf",
    width = 9,
    height = 4)
ggarrange(richness.plot,
          E.plot)
dev.off()
