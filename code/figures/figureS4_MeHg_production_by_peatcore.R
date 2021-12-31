#### code/figures/figureS3_MeHg_production_by_porewater.R ####
# Benjamin D. Peterson


#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")


#### Read in incubation data ####
inc.Hg.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds")


#### Order by average RMP ####
inc.Hg.data.order <- inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(RMP_core = mean(RMP_core)) %>%
  arrange(desc(RMP_core)) %>%
  select(coreID) %>%
  mutate(coreID = as.character(coreID)) %>%
  unlist()
inc.Hg.data <- inc.Hg.data %>%
  mutate(coreID = fct_relevel(coreID,
                                inc.Hg.data.order))
color.vector <- color.vector[inc.Hg.data.order]

#### Plot MeHg production, faceted by sediment core, ordered by RMP ####
inc.data.plot.by.core <- inc.Hg.data %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent) * 100,
            meth_spike_per_sd = sd(SMHG_201_percent),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  mutate(matrixID = paste("Pore water source: ", matrixID, sep = ""),
         matrixID = fct_relevel(matrixID, paste("Pore water source: ", MG.order, sep = ""))) %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = coreID)) +
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  facet_wrap(~matrixID, nrow = 2) +
  scale_fill_manual(values = color.vector) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black",
                                   size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())


#### Save out plot ####
pdf("results/figures/S4_MeHg_production_ordered_by_RMP_of_cores.pdf",
    width = 12,
    height = 6)
inc.data.plot.by.core
dev.off()
