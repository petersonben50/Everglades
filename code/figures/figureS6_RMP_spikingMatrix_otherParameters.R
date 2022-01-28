#### code/figures/figureS5_RMP_spikingMatrix_otherParameters.R ####
# Benjamin D. Peterson


#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Read in data
RMP.sulfide.log <- readRDS("results/incubations/RMP_porewater_sulfide.rds")
RMP.DOC <- readRDS("results/incubations/RMP_porewater_doc.rds")
RMP.UV <- readRDS("results/incubations/RMP_porewater_uv.rds")


#### Set it up ####
figures <- ggarrange(RMP.DOC,
                     RMP.UV,
                     RMP.sulfide.log,
                     ncol = 1)


#### Save out figure ####
pdf("results/figures/S6_RMP_spikingMatrix_otherParameters.pdf",
    height = 11,
    width = 5)
figures
dev.off()
