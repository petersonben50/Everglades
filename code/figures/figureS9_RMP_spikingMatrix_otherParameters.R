#### code/figures/figureS9_RMP_spikingMatrix_otherParameters.R ####
# Benjamin D. Peterson


#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Read in data
RMP.sulfide.log <- readRDS("results/incubations/RMP_porewater_sulfide.rds") + theme(legend.position = "none")
RMP.sulfate.log <- readRDS("results/incubations/RMP_porewater_sulfate.rds") + theme(legend.position = "none")
RMP.DOC <- readRDS("results/incubations/RMP_porewater_doc.rds") + theme(legend.position = "none")
RMP.UV <- readRDS("results/incubations/RMP_porewater_uv.rds")


#### Set it up ####
figures <- ggarrange(RMP.DOC,
                     RMP.UV,
                     RMP.sulfide.log,
                     RMP.sulfate.log,
                     ncol = 2,
                     nrow = 2,
                     labels = c("a", "b", "c", "d"))


#### Save out figure ####
pdf("results/figures/S9_RMP_spikingMatrix_otherParameters.pdf",
    height = 6,
    width = 7.2)
figures
dev.off()
