#### code/geochemistry/cleaning_geochem_2018.R ####
# Benjamin D. Peterson

# This script cleans up the geochemical data that 
# Brett sent me corresponding to the 2018 metagenomes.


#### Gettin' ready ####

rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(readxl)
library(tidyverse)


#### Read in data ####

raw.geochem.data <- read_xlsx("dataRaw/geochem/Data and figures for Ben Peterson_20120518_with sulfide.xlsx",
                              skip = 1) %>%
  rename(site = Site,
         PW_DOC_mg_L = `DOC [mg/L C]...3`,
         PW_SUVA_254 = `SUVA, 254nm...4`,
         PW_sulfate_mg_L = `Sulfate (mg/L)...5`,
         PW_sulfide_ug_L = `Sulfide (ug/L)...6`,
         PW_FTHg_ng_L = `FTHg (ng/L)...7`,
         PW_FMHg_ng_L = `F-MeHg (ng/L)...8`,
         SW_DOC_mg_L = `DOC [mg/L C]...10`,
         SW_SUVA_254 = `SUVA, 254nm...11`,
         SW_sulfate_mg_L = `Sulfate (mg/L)...12`,
         SW_sulfide_ug_L = `Sulfide (ug/L)...13`,
         SW_FTHg_ng_L = `FTHg (ng/L)...14`,
         SW_FMHg_ng_L = `F-MeHg (ng/L)...15`) %>%
  select(site,
         PW_DOC_mg_L,
         PW_SUVA_254,
         PW_sulfate_mg_L,
         PW_sulfide_ug_L,
         PW_FTHg_ng_L,
         PW_FMHg_ng_L,
         SW_DOC_mg_L,
         SW_SUVA_254,
         SW_sulfate_mg_L,
         SW_sulfide_ug_L,
         SW_FTHg_ng_L,
         SW_FMHg_ng_L) %>%
  filter(!is.na(site))


#### Read in WW data ####

raw.geochem.data.WW <- read_xlsx("dataRaw/geochem/Data and figures for Ben Peterson_20190415_with wagon wheel.xlsx",
                                 skip = 1) %>%
  rename(site = Site,
         PW_DOC_mg_L = `DOC [mg/L C]...3`,
         PW_SUVA_254 = `SUVA, 254nm...4`,
         PW_sulfate_mg_L = `Sulfate (mg/L)...5`,
         PW_sulfide_ug_L = `Sulfide (mg/L)...6`,
         PW_FTHg_ng_L = `FTHg (ng/L)...7`,
         PW_FMHg_ng_L = `F-MeHg (ng/L)...8`,
         SW_DOC_mg_L = `DOC [mg/L C]...10`,
         SW_SUVA_254 = `SUVA, 254nm...11`,
         SW_sulfate_mg_L = `Sulfate (mg/L)...12`,
         SW_sulfide_ug_L = `Sulfide (mg/L)...13`,
         SW_FTHg_ng_L = `FTHg (ng/L)...14`,
         SW_FMHg_ng_L = `F-MeHg (ng/L)...15`) %>%
  select(site,
         PW_DOC_mg_L,
         PW_SUVA_254,
         PW_sulfate_mg_L,
         PW_sulfide_ug_L,
         PW_FTHg_ng_L,
         PW_FMHg_ng_L,
         SW_DOC_mg_L,
         SW_SUVA_254,
         SW_sulfate_mg_L,
         SW_sulfide_ug_L,
         SW_FTHg_ng_L,
         SW_FMHg_ng_L) %>%
  filter(site == "Big Cypress Wagon Wheel")


#### Clean WW data ####
# Rename it
raw.geochem.data.WW$site <- "WW"
# The sulfide measurements at this site were not given to me as a concentration, just a large negative number.
# We'll take em out.
raw.geochem.data.WW$SW_sulfide_ug_L <- NA
raw.geochem.data.WW$PW_sulfide_ug_L <- NA


#### Combine data ####

raw.geochem.data <- rbind(raw.geochem.data,
                          raw.geochem.data.WW)


#### Clean up metadata columns ####

geochem.data <- raw.geochem.data %>%
  gather(key = constituent,
         value = concentration,
         -1 ) %>%
  mutate(location = strsplit(constituent, "_") %>% sapply("[", 1),
         constituent = sapply(strsplit(constituent, "_"),
                              function(x) {
                                paste(x[-1], collapse = '_')
                              })) %>%
  filter(!(site %in% c("E01", "T01", "UC1")))



#### Replace sulfate value of <0.2 with 0 ####

geochem.data$concentration[which(geochem.data$concentration == '<0.2')] <- 0
geochem.data <- geochem.data %>%
  mutate(concentration = as.numeric(concentration),
         concentration = round(concentration, 3))

 
#### Make data table wide ####
geochem.data <- geochem.data %>%
  spread(key = constituent,
         value = concentration) %>%
  arrange(location)


#### Read out cleaned data ####

write.csv(geochem.data,
          "dataEdited/geochem/geochem_data_2018.csv",
          row.names = FALSE)
