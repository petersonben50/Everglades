#### code/assembly_analysis/hgcA/hgcA_mcrA_ratio.R ####
# Benjamin D. Peterson



#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")



##### Read in MG metadata ####
MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.renaming.vector <- MG.metadata$siteID
names(site.renaming.vector) <- MG.metadata$metagenomeID

medium.renaming.vector <- MG.metadata$medium
names(medium.renaming.vector) <- MG.metadata$metagenomeID




#### Read in methanogen hgcA data ####
methanogen.hgcA.data <- read.csv("dataEdited/assembly_analysis/hgcA/depth/hgcA_coverage_scgNormalization.csv",
                                 stringsAsFactors = FALSE) %>%
  full_join(read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                      sheet = "seq_class")) %>%
  arrange(classification) %>%
  filter(seqID %in% readLines("dataEdited/assembly_analysis/hgcA/hgcA_true.txt"))%>%
  filter(clusterName == "Methanogen")

methanogen.hgcA.data.df <- methanogen.hgcA.data %>%
  select(grep("KMBP005", names(methanogen.hgcA.data))) %>%
  gather(key = metagenomeID,
         value = hgcA.abundance) %>%
  group_by(metagenomeID) %>%
  summarise(hgcA.abundance = sum(hgcA.abundance))




##### Read in methanogen info ####

# Read in mcrA list
mcrA.scaffold.list <- readLines("dataEdited/assembly_analysis/metabolicProteins/methanogenesis/mcrA/mcrA_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)

methanogen.mcrA.data <- read.csv("dataEdited/assembly_analysis/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                                 stringsAsFactors = FALSE) %>%
  filter(scaffoldID %in% mcrA.scaffold.list) 

methanogen.mcrA.data.df <- methanogen.mcrA.data %>%
  select(grep("KMBP005", names(methanogen.mcrA.data))) %>%
  gather(key = metagenomeID,
       value = mcrA.abundance) %>%
  group_by(metagenomeID) %>%
  summarise(mcrA.abundance = sum(mcrA.abundance))



#### Combine methanogen data ####
methanogen.data <- full_join(methanogen.hgcA.data.df,
                             methanogen.mcrA.data.df) %>%
  mutate(siteID = site.renaming.vector[metagenomeID],
         medium = medium.renaming.vector[metagenomeID]) %>%
  mutate(hgcA_mcrA_percent = hgcA.abundance / mcrA.abundance * 100) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  arrange(hgcA.abundance)
methanogen.data
