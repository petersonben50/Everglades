#### code/2019_analysis_assembly/non_respiratory/rnf_abundance.R ####
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





#### Read in marker lists ####

list.of.marker.lists <- list.files(path = "dataEdited/2019_analysis_assembly/metabolicProteins/non_respiration/rnf/",
                                   pattern = "_derep_list.txt")
for (marker.of.interest in 1:length(list.of.marker.lists)) {
  scaffolds.of.interest = readLines(paste("dataEdited/2019_analysis_assembly/metabolicProteins/non_respiration/rnf/",
                                          list.of.marker.lists[marker.of.interest],
                                          sep = "")) %>%
    strsplit("_[1-9]+") %>%
    sapply("[", 1)
  geneNameOfInterest = gsub("_derep_list.txt", "", list.of.marker.lists[marker.of.interest])
  marker.df.temp <- data.frame(scaffoldID = scaffolds.of.interest,
                               geneName = rep(geneNameOfInterest,
                                              length(scaffolds.of.interest)))
  
  if (marker.of.interest == 1) {
    marker.df <- marker.df.temp
  } else {
    marker.df <- rbind(marker.df,
                       marker.df.temp)
  }
  rm(marker.df.temp)
}

marker.depth <- left_join(marker.df,
                          depth.data)
unique(marker.depth$geneName)



#### Generate plot ####

color.vector.rnf <- cb.translator[1:length(unique(marker.depth$geneName))]
names(color.vector.rnf) <- unique(marker.depth$geneName)
rnf.plot <- plot.scaffold.coverage(marker.depth.df = marker.depth,
                                   genesOfInterest = unique(marker.depth$geneName),
                                   medium.of.interest = "sediment",
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.rnf,
                                   MG.order.vector = MG.order,
                                   titleToUse = paste("Rnf complex subunits"))
pdf("results/2019_analysis_assembly/non_respiration/rnf_abundance.pdf",
    height = 3,
    width = 6)
rnf.plot
dev.off()
