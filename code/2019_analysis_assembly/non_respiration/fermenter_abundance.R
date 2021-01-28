#### code/2019_analysis_assembly/non_respiration/fermenter_abundance.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(ggplot2)
library(patchwork)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/2019_analysis_assembly/metabolic_protein_plots.R")


#### Set site order for plotting ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")


##### Read in MG metadata ####
MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.renaming.vector <- MG.metadata$siteID
names(site.renaming.vector) <- MG.metadata$metagenomeID
medium.renaming.vector <- MG.metadata$medium
names(medium.renaming.vector) <- MG.metadata$metagenomeID






##### Read in depth info ####

depth.data <- read.csv("dataEdited/2019_analysis_assembly/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE) %>%
  gather(key = sampleID,
         value = coverage,
         -1) %>%
  mutate(siteID = site.renaming.vector[sampleID]) %>%
  mutate(medium = medium.renaming.vector[sampleID])


#### Plot out PFOR (this is really all I got) ####
pfor.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/non_respiration/pfor_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
pfor.tree.viz <- plot.scaffold.coverage.for.one.gene(depth.data.of.interest = depth.data,
                                                     scaffold.list = pfor.scaffolds,
                                                     geneName = "pfor",
                                                     medium.of.interest = "sediment",
                                                     use.points = FALSE)

pdf("results/2019_analysis_assembly/non_respiration/pfor_abundance.pdf",
    height = 3,
    width = 6)
pfor.tree.viz
dev.off()
