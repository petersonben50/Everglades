#### code/incubations/2019_inc_amb.R ####
# Benjamin D. Peterson

#### Get cleaned up ####
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
amb.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis of core experiments_v2.xlsx",
                         sheet = "incubation_results_cleaned",) %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) 


#### Generate plot of ambient MeHg ####
pdf("results/2019_incubations/ambient_MeHg.pdf",
    width = 6,
    height = 3)
amb.Hg.data %>%
  group_by(coreID) %>%
  summarise(MeHg_amb_mean = mean(MeHg_amb),
            MeHg_amb_sd = sd(MeHg_amb),
            MeHg_amb_count = n(),
            MeHg_amb_se = MeHg_amb_sd / sqrt(MeHg_amb_count)) %>%
  ggplot(aes(x = coreID,
             y = MeHg_amb_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = MeHg_amb_mean - MeHg_amb_se,
                    ymax = MeHg_amb_mean + MeHg_amb_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector.core) +
  labs(title = "Ambient MeHg in cores",
       y = "MeHg (ng/L)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
dev.off()



#### Save out ambient MeHg data ####
amb.Hg.data %>%
  saveRDS("dataEdited/2019_incubations/ambient_MeHg.rds")


