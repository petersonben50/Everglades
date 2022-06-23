#### code/figures/figure2_RMP_porewater_suva.R ####
# Benjamin D. Peterson

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(lme4)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")



#### Read in plots ####
RMP.SUVA <- readRDS("results/incubations/RMP_porewater_SUVA.rds")


#### Print out graphs ####
pdf("results/figures/2_RMP_porewater.pdf",
    height = 6,
    width = 7)
RMP.SUVA
dev.off()
