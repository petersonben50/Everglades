#### code/figures/figure5_RMP_ambient_MeHg.R ####
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
RMP.amb.MeHg.prod <- readRDS("dataEdited/incubations/RMP_ambient_MeHg_production_comparison.rds")
sed.MeHg.vs.RMP.core <- readRDS("results/incubations/ambient_MeHg_vs_RMP_core.rds")
sed.MeHg.vs.RMP.matrix <- readRDS("results/incubations/ambient_MeHg_vs_RMP_matrix.rds")




#### 3D scatterplot of RMP and native MeHg production ####
pdf("results/figures/5/3D_plot_RMP_native_methylation.pdf",
    height = 4.5,
    width = 6)
scatterplot3d(data.frame(RMP.amb.MeHg.prod$RMP_core_mean,
                         RMP.amb.MeHg.prod$RMP_porewater_mean,
                         RMP.amb.MeHg.prod$SMHG_201_mean),
              angle = 40,
              pch = 16,
              cex.symbols = 1.5,
              color = color.vector[RMP.amb.MeHg.prod$siteID],
              type = "h",
              lab = c(4, 4),
              lab.z = 4,
              xlim = c(0, 80),
              ylim = c(0, 80),
              zlim = c(0, 3),
              xlab = "Core RMP",
              ylab = "Porewater RMP",
              zlab = "MeHg production\nunder native conditions")
legend(x = 0.6, y = 4.7,
       legend = RMP.amb.MeHg.prod$siteID,
       pch = 16,
       col = color.vector[RMP.amb.MeHg.prod$siteID])

dev.off()


#### Plots of RMP values vs. ambient levels ####
pdf("results/figures/5/RMP_vs_ambient_MeHg.pdf",
    height = 4,
    width = 10)

ggarrange(sed.MeHg.vs.RMP.matrix + theme(legend.position = "none"),
          sed.MeHg.vs.RMP.core + theme(legend.position = "none"),
          ncol = 2)
dev.off()
