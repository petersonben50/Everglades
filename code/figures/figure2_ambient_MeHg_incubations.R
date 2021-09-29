#### code/incubations/figure2_ambient_MeHg_incubations.R ####
# Benjamin D. Peterson


#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(lme4)
library(patchwork)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")



#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                         sheet = "cleaned_incubation_data",) %>%
  mutate(coreID = fct_relevel(coreID, MG.order),
         matrixID = fct_relevel(matrixID, PW.order))



#### Read in porewater data, order it, and generate percent MeHg ####
porewater.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                            sheet = "porewater_data") %>%
  filter(siteID %in% MG.order) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  mutate(MeHg_porewater = as.numeric(FMHG),
         HgT_porewater = as.numeric(FTHG),
         percent_MeHg_porewater = MeHg_porewater / HgT_porewater * 100) %>%
  select(siteID, MeHg_porewater, HgT_porewater, percent_MeHg_porewater)



#### Calculate ambient MeHg in sediment from incubation cores ####
amb.Hg.sediment <- inc.Hg.data %>%
  group_by(coreID) %>%
  mutate(SMHG_201_percent = SMHG_201_percent * 100) %>%
  summarise(SMHG_amb_mean = mean(SMHG_amb),
            SMHG_amb_sd = sd(SMHG_amb),
            SMHG_amb_count = n(),
            SMHG_amb_se = SMHG_amb_sd / sqrt(SMHG_amb_count),
            STHG_amb_mean = mean(STHG_amb),
            STHG_amb_sd = sd(STHG_amb),
            STHG_amb_count = n(),
            STHG_amb_se = STHG_amb_sd / sqrt(STHG_amb_count),
            SMHG_amb_percent_mean = mean(SMHG_amb_percent),
            SMHG_amb_percent_sd = sd(SMHG_amb_percent),
            SMHG_amb_percent_count = n(),
            SMHG_amb_percent_se = SMHG_amb_percent_sd / sqrt(SMHG_amb_percent_count)) %>%
  select(coreID, SMHG_amb_mean, SMHG_amb_se,
         STHG_amb_mean, STHG_amb_se,
         SMHG_amb_percent_mean, SMHG_amb_percent_se)
amb.Hg.sediment



#### Generate plot of ambient MeHg in soil ####
MeHg.soil.plot <- amb.Hg.sediment %>%
  ggplot(aes(x = coreID,
             y = SMHG_amb_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = SMHG_amb_mean - SMHG_amb_se,
                    ymax = SMHG_amb_mean + SMHG_amb_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(y = "Sediment MeHg (ng/g)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 6),
        legend.position = "none",
        plot.margin = margin(b = 10,
                             l = 25,
                             t = 10,
                             r = 10))
MeHg.soil.plot


#### Porewater MeHg concentration plot ####
MeHg.porewater.plot <- porewater.data %>%
  ggplot(aes(x = siteID,
             y = MeHg_porewater,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector) +
  labs(y = "Porewater MeHg (ng/L)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(b = 10,
                             l = 13,
                             t = 10,
                             r = 10))
MeHg.porewater.plot



#### Plot MeHg production across cores, points linked by spiking matrix ####
MeHg.production.data <- inc.Hg.data %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent) * 100,
            meth_spike_per_sd = sd(SMHG_201_percent),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  mutate(native_PW = (as.character(coreID) == as.character(matrixID)) * 1)

MeHg.production.plot <- MeHg.production.data %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             group = matrixID,
             col = matrixID)) +
  geom_point(aes(x = coreID,
                 y = meth_spike_per_mean,
                 stroke = 1.5),
             data = MeHg.production.data %>%
               filter(native_PW == 1),
             shape = 4,
             size = 3,
             color = "black") +
  geom_point() +
  geom_line() +
  scale_color_manual(name = "Spiking matrix",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank(),
        legend.position = c(0.15, 0.68))
MeHg.production.plot


#### Set up ordering of plot ####
figure <- ggarrange(ggarrange(MeHg.soil.plot,
                              MeHg.porewater.plot,
                              nrow = 2,
                              labels = c("A.", "B."),
                              label.x = -0.02,
                              label.y = 1.02),
                    MeHg.production.plot,
                    labels = c("", "C."),
                    ncol = 2,
                    widths = c(1,2),
                    label.x = -0.01,
                    label.y = 1.01)
figure



#### Save out plot ####
pdf("results/incubations/figure_2.pdf",
    width = 9,
    height = 5)
figure
dev.off()
