#### code/incubations/native_MeHg_production_function_of_RMP.R ####
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
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")



#### Read in incubation data ####
Hg.inc.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds")
color.vector <- color.vector[unique(Hg.inc.data$coreID)]


#### Calculate core RMP mean and SE ####
coreID.data <- Hg.inc.data %>%
  group_by(coreID) %>%
  summarise(RMP_core_mean = mean(RMP_core),
            RMP_core_sd = sd(RMP_core),
            RMP_core_count = n(),
            RMP_core_se = RMP_core_sd / sqrt(RMP_core_count)) %>%
  ungroup() %>%
  rename(siteID = coreID) %>%
  select(siteID, RMP_core_mean, RMP_core_sd, RMP_core_se)


#### Calculate porewater RMP mean and SE ####
matrixID.data <- Hg.inc.data %>%
  group_by(matrixID) %>%
  summarise(RMP_porewater_mean = mean(RMP_porewater),
            RMP_porewater_sd = sd(RMP_porewater),
            RMP_porewater_count = n(),
            RMP_porewater_se = RMP_porewater_sd / sqrt(RMP_porewater_count)) %>%
  ungroup() %>%
  filter(matrixID %in% MG.order) %>%
  rename(siteID = matrixID) %>%
  select(siteID, RMP_porewater_mean, RMP_porewater_sd, RMP_porewater_se)


#### Calculate ambient MeHg mean and SE ####
ambient.data <- Hg.inc.data %>%
  group_by(coreID) %>%
  summarise(amb_MeHg_mean = mean(SMHG_amb_percent),
            amb_MeHg_sd = sd(SMHG_amb_percent),
            amb_MeHg_count = n(),
            amb_MeHg_se = amb_MeHg_sd / sqrt(amb_MeHg_count)) %>%
  ungroup() %>%
  rename(siteID = coreID) %>%
  select(siteID, amb_MeHg_mean, amb_MeHg_sd, amb_MeHg_se)



#### Read in porewater data, order it, and generate percent MeHg ####
porewater.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                            sheet = "porewater_data") %>%
  filter(siteID %in% MG.order) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  mutate(MeHg_porewater = as.numeric(FMHG),
         HgT_porewater = as.numeric(FTHG),
         percent_MeHg_porewater = MeHg_porewater / HgT_porewater * 100) %>%
  select(siteID, MeHg_porewater, HgT_porewater, percent_MeHg_porewater)




#### Combine data ####
all.data <- coreID.data %>%
  full_join(matrixID.data) %>%
  full_join(ambient.data) %>%
  full_join(porewater.data)



#### Plot: Ambient sediment MeHg vs. RMP values ####
sed.MeHg.vs.RMP.core <- all.data %>%
  ggplot(aes(x = RMP_core_mean,
             y = amb_MeHg_mean,
             color = siteID)) +
  geom_point(size = 4) +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  theme_bw() +
  xlim(c(0, 80)) +
  ylab("Ambient MeHg (ng/L)") +
  xlab("Peat core RMP") +
  scale_y_continuous(trans = 'log10',
                     limits = c(0.001, 0.1)) +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = c(0.2, 0.6))

sed.MeHg.vs.RMP.matrix <- all.data %>%
  ggplot(aes(x = RMP_porewater_mean,
             y = amb_MeHg_mean,
             color = siteID)) +
  geom_point(size = 4) +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  theme_bw() +
  xlim(c(0, 80)) +
  ylab("Ambient MeHg (ng/g)") +
  xlab("Pore water RMP") +
  scale_y_continuous(trans = 'log10',
                     limits = c(0.001, 0.1)) +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = c(0.2, 0.6))
# Plot them together
ggarrange(sed.MeHg.vs.RMP.core,
          sed.MeHg.vs.RMP.matrix + theme(legend.position = "none"),
          labels = c("A.", "B."),
          ncol = 2)

# Save out figures for later use in figures
saveRDS(sed.MeHg.vs.RMP.core,
        "results/incubations/ambient_MeHg_vs_RMP_core.rds")
saveRDS(sed.MeHg.vs.RMP.matrix,
        "results/incubations/ambient_MeHg_vs_RMP_matrix.rds")


#### Stats: Ambient sediment MeHg vs. RMP values ####
ambient.sed.vs.RMPpeat.model <- lm(log(amb_MeHg_mean, 10) ~ RMP_core_mean,
                                   data = all.data)
summary(ambient.sed.vs.RMPpeat.model)

ambient.sed.vs.RMPmatrix.model <- lm(amb_MeHg_mean ~ RMP_porewater_mean,
                                     data = all.data)
summary(ambient.sed.vs.RMPmatrix.model)



#### Plot: Ambient sediment MeHg vs. RMPpeat with linear fit ####
sed.MeHg.vs.RMP.core.model <- all.data %>%
  ggplot(aes(x = RMP_core_mean,
             y = amb_MeHg_mean,
             color = siteID)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  geom_point(size = 4) +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  theme_bw() +
  xlim(c(0, 100)) +
  ylab("Ambient MeHg (ng/L)") +
  xlab("Peat core RMP") +
  scale_y_continuous(trans = 'log10',
                     limits = c(0.001, 0.1)) +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = c(0.2, 0.6))

pdf("results/figures/S13_ambientMeHg_vs_RMP.pdf",
    width = 8,
    height = 4)
ggarrange(sed.MeHg.vs.RMP.core.model,
          sed.MeHg.vs.RMP.matrix + theme(legend.position = "none"),
          labels = c("A.", "B."),
          ncol = 2)
dev.off()

#### Plot: Ambient pore water MeHg against the RMP values ####
PW.MeHg.vs.RMP.core <- all.data %>%
  ggplot(aes(x = RMP_core_mean,
             y = log(percent_MeHg_porewater),
             color = siteID)) +
  geom_point(size = 3) +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  theme_bw() +
  xlim(c(0, 80)) +
  ylab("Ambient porewater MeHg log(ng/L)") +
  xlab("Peat core RMP") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.2, 0.6))

PW.MeHg.vs.RMP.matrix <- all.data %>%
  ggplot(aes(x = RMP_porewater_mean,
             y = log(percent_MeHg_porewater),
             color = siteID)) +
  geom_point(size = 3) +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  theme_bw() +
  xlim(c(0, 80)) +
  ylab("Ambient porewater MeHg log(ng/L)") +
  xlab("Pore water RMP") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.2, 0.6))

ggarrange(PW.MeHg.vs.RMP.core,
          PW.MeHg.vs.RMP.matrix + theme(legend.position = "none"),
          ncol = 2)
