#### code/figures/figureS12_ambientMeHg_vs_RMP.R ####
# Benjamin D. Peterson



#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(scatterplot3d)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")




#### Read in data ####
sed.MeHg.vs.RMP.core <- readRDS("results/incubations/ambient_MeHg_vs_RMP_core.rds")
sed.MeHg.vs.RMP.matrix <- readRDS("results/incubations/ambient_MeHg_vs_RMP_matrix.rds")




#### Plots of RMP values vs. ambient levels ####
pdf("results/figures/S12_ambientMeHg_vs_RMP.pdf",
    height = 4,
    width = 10)

ggarrange(sed.MeHg.vs.RMP.matrix + theme(legend.position = "none"),
          sed.MeHg.vs.RMP.core + theme(legend.position = "none"),
          labels = c("A.", "B."),
          ncol = 2)
dev.off()
