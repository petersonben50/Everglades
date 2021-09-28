#### code/incubations/figure3_RMP_covariates.R ####
# Benjamin D. Peterson

# Obsidian notes here: Incubations - statistical analysis porewaters
# Results: results/RMP_porewater.pptx

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")



#### Read in plots ####
RMP.hgcA <- readRDS("results/incubations/RMP_microbes_hgcA.rds")
RMP.SUVA <- readRDS("results/incubations/RMP_porewater_SUVA.rds")



#### Set up ordering of plot ####
figure <- ggarrange(RMP.SUVA,
                    RMP.hgcA,
                    labels = c("A.", "B."),
                    ncol = 2)



#### Print out graphs ####
pdf("results/incubations/figure_3_RMP.pdf",
    height = 5.5,
    width = 12)
figure
dev.off()

