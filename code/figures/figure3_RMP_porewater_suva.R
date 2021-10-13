#### code/figures/figure3_RMP_porewater_suva.R ####
# Benjamin D. Peterson

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
RMP.SUVA <- readRDS("results/incubations/RMP_porewater_SUVA.rds")


#### Print out graphs ####
pdf("results/figures/3_RMP_porewater.pdf",
    height = 6,
    width = 7)
RMP.SUVA
dev.off()
