#### code/geochemistry/2019_geochem.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(patchwork)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in incubation data ####
inc.Hg.data <- read_xlsx("dataEdited/geochem/geochem_data_2019.xlsx",
                         sheet = "incubations")


#### Order of samples to be plotted ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")


#### Read in PW geochem data ####
all.data <- read_xlsx("dataEdited/geochem/geochem_data_2019.xlsx",
                      sheet = "PW_SW_geochem") %>%
  mutate(siteID = fct_relevel(siteID, MG.order))


#### Prep color vector ####
color.vector.core <- cb.translator[1:6]
names(color.vector.core) <- MG.order


#### Generate plots for figure 1 ####
sulfate.SW.plot <- all.data %>%
  filter(medium == "SW") %>%
  ggplot(aes(x = siteID,
             y = Sulfate_mg.L / 96.06 * 1000, # Divide by molar mass and multiple by 1000 to get µM
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector.core) +
  ylim(c(0, 450)) +
  labs(y = "Sulfate (µM)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")
sulfate.PW.plot <- all.data %>%
  filter(medium == "PW") %>%
  ggplot(aes(x = siteID,
             y = Sulfate_mg.L / 96.06 * 1000, # Divide by molar mass and multiple by 1000 to get µM
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector.core) +
  ylim(c(0, 190)) +
  labs(y = "Sulfate (µM)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")
sulfide.PW.plot <- all.data %>%
  filter(medium == "PW") %>%
  ggplot(aes(x = siteID,
             y = Sulfide_µg.L / 34.08, # Divide by molar mass and multiple by 1000 to get µM
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector.core) +
  ylim(c(0, 110)) +
  labs(y = "Sulfide (µM)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")


#### Save out plots for figure 1 ####
pdf("results/geochemistry/2019_S_manuscript.pdf",
    width = 4,
    height = 5)
sulfate.SW.plot / sulfate.PW.plot / sulfide.PW.plot
dev.off()


#### Generate plots for Supplementary Figure 1 ####
sulfide.SW.plot <- all.data %>%
  filter(medium == "SW") %>%
  ggplot(aes(x = siteID,
             y = Sulfide_µg.L / 34.08, # Divide by molar mass and multiple by 1000 to get µM
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector.core) +
  ylim(c(0, 30)) +
  labs(y = "Sulfide (µM)",
       title = "A.") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none") +
  geom_text(x = 5, y = 3,
            label = "n.d.",
            fontface = "plain") +
  geom_text(x = 6, y = 3,
            label = "n.d.",
            fontface = "plain")

DOC.SW.plot <- all.data %>%
  filter(medium == "SW") %>%
  ggplot(aes(x = siteID,
             y = as.numeric(DOC_mg.L),
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector.core) +
  labs(y = "DOC (mg/L)",
       title = "B.") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")

suva.SW.plot <- all.data %>%
  filter(medium == "SW") %>%
  ggplot(aes(x = siteID,
             y = as.numeric(SUVA_254_L.mg.m),
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector.core) +
  ylim(c(0, 4)) +
  labs(y = "SUVA (L/mg/m)",
       title = "C.") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_blank()) +
  theme(legend.position = "none")


#### Save out plots for Supplementary Figure 1 ####
pdf("results/geochemistry/2019_BGC_manuscriptSupplement.pdf",
    width = 4,
    height = 5)
sulfide.SW.plot / DOC.SW.plot / suva.SW.plot
dev.off()
