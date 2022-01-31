#### code/figures/figureS10_tax_rel_coverage.R ####
# Benjamin D. Peterson


#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Read in ggplot of graph ####
tax.group.coverage.relative.plot <- readRDS("results/metagenomes/assembly/hgcA/tax_group_rel_coverage.rds")



#### Save out figure ####
pdf("results/figures/S10_tax_rel_coverage.pdf",
    width = 7.5,
    height = 4)
tax.group.coverage.relative.plot
dev.off()
