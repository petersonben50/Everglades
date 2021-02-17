rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Generate vector of correct order of samples along sulfate gradient ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")


#### Prep color vector ####
color.vector.core <- cb.translator[1:6]
names(color.vector.core) <- MG.order


#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis of core experiments_v2.xlsx",
                         sheet = "incubation_results_cleaned") %>%
  filter(coreID == matrixID) %>%
  mutate(coreID = fct_relevel(coreID, MG.order))


#### Generate plot #### 
pdf("results/2019_incubations/MeHg_production_native.pdf",
    width = 6,
    height = 3)
inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(meth_spike_per_mean = mean(MeHg_fract_spike * 100),
            meth_spike_per_sd = sd(MeHg_fract_spike * 100),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count)) %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector.core) +
  labs(y = "Methylated spike (%)",
       title = "Hg methylation using ambient porewater") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
dev.off()


#### Save out methylation data under native conditions ####
inc.Hg.data %>%
  saveRDS("dataEdited/2019_incubations/native_methylation.rds")
