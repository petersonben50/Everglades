#### code/geochemistry/clean_geochem_PW.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")


#### Read in PW geochem data ####
all.PW.data <- read_excel("dataRaw/geochem/2019/Everglades_WCA_data_release_v1.4_FINAL.xlsm",
                          sheet = "Table_3_Water") %>%
  filter(`Sample Collection Code` == 'b',
         year(`Sample Collection Date (mm/dd/yy h:mm)`) == 2019) %>%
  mutate(siteID = gsub(" WI TRANSECT [1:2]-", "-",
                       `Site Name`)) %>%
  rename(sampleDate = `Sample Collection Date (mm/dd/yy h:mm)`,
         FTHG = `f.THg (ng/L)`,
         FMHG = `f.MeHg (ng/L)`,
         PTHG = `p.THg (ng/L)`,
         PMHG = `p.MeHg (ng/L)`,
         DOC = `DOC (mgC/L)`,
         UVabs = `UV Abs 254nm`,
         SUVA = `SUVA254nm (L/mgC m)`,
         Cl = `Chloride (mg/L)`,
         Br = `Bromide (mg/L)`,
         sulfate = `Sulfate (mg/L)`,
         Li = `Lithium (mg/L)`,
         Na = `Sodium (mg/L)`,
         K = `Potassium (mg/L)`,
         Mg = `Magnesium (mg/L)`,
         Ca = `Calcium (mg/L)`,
         sulfide = `Sulfide (µg/L)`,
         ORP = `ORP (mV)`,
         conductivity = `Cond (µS/cm)`,
         DO = `DO (mg/L)`) %>%
  select(siteID, sampleDate,
         FTHG, FMHG, PTHG, PMHG,
         DOC, UVabs, SUVA, sulfate, sulfide,
         Cl, Br, Li, Na, K, Mg, Ca,
         ORP, conductivity, DO, pH) %>%
  mutate(sulfide = as.numeric(gsub("--", 0, sulfide)) / 1000) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:2)) %>%
  mutate(concentration = as.numeric(gsub("--", 0, concentration)))


#### Add in the modified data ####
all.PW.data[((all.PW.data$siteID == "2A-N") &
               (all.PW.data$constituent == "DOC")),
            "concentration"] <- 49

all.PW.data[((all.PW.data$siteID == "2A-N") &
               (all.PW.data$constituent == "SUVA")),
            "concentration"] <- 3.66

all.PW.data[((all.PW.data$siteID == "2A-N") &
               (all.PW.data$constituent == "UVabs")),
            "concentration"] <- 1.79




#### Save out data ####

write.csv(x = all.PW.data,
          "dataEdited/geochem/clean_porewater_2019_geochem.csv",
          row.names = FALSE)
