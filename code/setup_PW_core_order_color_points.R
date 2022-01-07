#### code/setup_PW_core_order_color_points.R ####
# Benjamin D. Peterson

# Additional notes in the corresponding Obsidian document:
# Statistical analysis of MeHg production - re-analyzed for paper


#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Generate vector of correct order of samples along sulfate gradient ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")
PW.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8", "F1", "CYSTEINE", "CONTROL")



#### Prep color and point vector ####
color.vector <- c(cb.translator[c("bluishgreen", "skyblue", "vermillion", "reddishpurple", "orange", "yellow", "black")],
                  "#7f7f7f", "#bfbfbf")
names(color.vector) <- PW.order
point.vector <- c(9, 18:15, 25, 8, 3, 14)
names(point.vector) <- PW.order



#### Add color info to GPS coordinates ####
# Read in GPS coordinates
gps.coord <- read_xlsx("dataEdited/gisData/EG_sampling_sites_GPS_2019.xlsx") %>%
  mutate(Site = Site %>%
           str_replace("WCA-", "")) %>%
  filter(Site %in% MG.order) %>%
  mutate(colorToUse = color.vector[Site])
write.csv(gps.coord,
          "dataEdited/gisData/EG_sampling_sites_GPS_2019.csv",
          row.names = FALSE)
