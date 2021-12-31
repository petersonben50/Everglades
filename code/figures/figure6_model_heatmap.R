#### code/figures/figure6_model_heatmap.R ####
# Benjamin D. Peterson

# Obsidian notes here: Incubations - statistical analysis porewaters
# Results: results/RMP_porewater.pptx

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(lme4)
library(readxl)
library(scatterplot3d)
library(tidyverse)
library(viridis)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")


colors.to.use <- colorRampPalette(colors = c(cb.translator["black"],
                                             cb.translator["bluishgreen"],
                                             cb.translator["vermillion"],
                                             cb.translator["yellow"]))

#### Prep data ####
model.data <- data.frame(x = rep(seq(0, 1, by = 0.01), 101) %>% sort(),
                         y = rep(seq(0, 1, by = 0.01), 101))
model.data$z <- (model.data$x) * (model.data$y)
# model.data$z <- (model.data$x)^2 * (model.data$y)^2

x <- seq(min(model.data$x),
         max(model.data$x),
         .01)
y <- seq(min(model.data$y),
         max(model.data$y),
         .01)

#### Read out pdf of heatmap ####

image(x = x,
      y = y,
      z = matrix(model.data$z,
                 ncol = 101),
      col = magma(100),
      breaks = seq(0, 1,
                   by = 0.01))
legend(col = magma(100))


pdf("results/figures/6/heatmap_MeHg_production.pdf",
    height = 6,
    width = 6)
filled.contour(x = x,
               y = y,
               z = matrix(model.data$z,
                          ncol = 101),
               xlim = c(0, 1),
               ylim = c(0, 1),
               zlim = c(0, 1),
               xaxt = "n",
               yaxt = "n",
               nlevels = 100,
               asp = 1,
               color.palette = magma,
               key.title = {par(cex.main=1);title(main="Relative\nMeHg Production")},
               key.axes = FALSE,
               frame.plot = FALSE)
dev.off()
