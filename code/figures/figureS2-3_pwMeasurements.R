#### code/geochemistry/geochem_PW.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(patchwork)
library(ggpubr)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")



#### Read in PW geochem data ####
all.data <- read.csv("dataEdited/geochem/clean_porewater_2019_geochem.csv") %>%
  mutate(siteID = fct_relevel(siteID,
                              MG.order))


#### Plot-maker function ####
geochem.plot.maker <- function(dataset = all.data,
                               constituent.of.interest,
                               y.label.to.use,
                               ylim.to.use = NULL,
                               color.vector.to.use = color.vector) {
  plot.plot.plot <- all.data %>%
    filter(constituent == constituent.of.interest) %>%
    ggplot(aes(x = siteID,
               y = concentration,
               width = 0.8,
               fill = siteID)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(values = color.vector.to.use) +
    labs(y = y.label.to.use) +
    theme_bw() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_blank()) +
    theme(legend.position = "none")

  if (!is.null(ylim.to.use)) {
    plot.plot.plot <- plot.plot.plot +
      ylim(ylim.to.use)
  }
  plot.plot.plot
}



#### Generate plots for figure 1 ####
sulfate.plot <- geochem.plot.maker(dataset = all.data,
                                   constituent.of.interest = "sulfate",
                                   y.label.to.use = "Sulfate (mg/L)")

sulfide.plot <- geochem.plot.maker(dataset = all.data,
                                   constituent.of.interest = "sulfide",
                                   y.label.to.use = "Sulfide (mg/L)")

DOC.plot <- geochem.plot.maker(dataset = all.data,
                               constituent.of.interest = "DOC",
                               y.label.to.use = "DOC (mgC/L)")

SUVA.plot <- geochem.plot.maker(dataset = all.data,
                                constituent.of.interest = "SUVA",
                                y.label.to.use = "SUVA (L/mgC/m)")

UVabs.plot <- geochem.plot.maker(dataset = all.data,
                                 constituent.of.interest = "UVabs",
                                 y.label.to.use = "UVabs (cm-1)")

Cl.plot <- geochem.plot.maker(dataset = all.data,
                              constituent.of.interest = "Cl",
                              y.label.to.use = "Cl (mg/L)")

#### Arrange plots ####
pdf("results/figures/geochem_plots.pdf",
    height = 7.5,
    width = 7)
ggarrange(sulfate.plot, sulfide.plot,
          DOC.plot, SUVA.plot,
          UVabs.plot, Cl.plot,
          labels = "AUTO",
          nrow = 3,
          ncol = 2)
dev.off()


#### Flow-through cell measurements ####
ORP.plot <- geochem.plot.maker(dataset = all.data,
                               constituent.of.interest = "ORP",
                               y.label.to.use = "ORP (mV)",
                               ylim.to.use = c(-180, 0))

DO.plot <- geochem.plot.maker(dataset = all.data,
                              constituent.of.interest = "DO",
                              y.label.to.use = "DO (mg/L)",
                              ylim.to.use = c(0, 1))

pH.plot <- geochem.plot.maker(dataset = all.data,
                              constituent.of.interest = "pH",
                              y.label.to.use = "pH",
                              ylim.to.use = c(0,10))

cond.plot <- geochem.plot.maker(dataset = all.data,
                                constituent.of.interest = "conductivity",
                                y.label.to.use = "Conductivity (µS/cm)",
                                ylim.to.use = c(0, 1200))


#### Arrange plots ####
pdf("results/figures/flowthrough_plots.pdf",
    height = 5,
    width = 7)
ggarrange(cond.plot, pH.plot,
          DO.plot, ORP.plot,
          labels = c("a", "b", "c", "d"),
          nrow = 2,
          ncol = 2)
dev.off()
