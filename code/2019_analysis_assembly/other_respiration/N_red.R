#### code/2019_analysis_assembly/other_respiration/n_reduction.R ####
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



#### Read in gene lists ####
list.of.marker.lists <- list.files(path = "dataEdited/2019_analysis_assembly/metabolicProteins/other_respiration/N/",
                                   pattern = "_derep_list.txt")
for (marker.of.interest in 1:length(list.of.marker.lists)) {
  
  scaffolds.of.interest = readLines(paste("dataEdited/2019_analysis_assembly/metabolicProteins/other_respiration/N/",
                                          list.of.marker.lists[marker.of.interest],
                                          sep = "")) %>%
    strsplit("_[1-9]+") %>%
    sapply("[", 1)
  gene_name <- gsub("_derep_list.txt", "", list.of.marker.lists[marker.of.interest])
  marker.df.temp <- data.frame(scaffoldID = scaffolds.of.interest,
                               geneName = rep(gene_name,
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



#### cox plots ####
color.vector.nrf <- c(cb.translator["vermillion"], cb.translator["yellow"])
names(color.vector.nrf) <- c("nrfA", "nrfH")
nrf.plot <- plot.scaffold.coverage(genesOfInterest = c("nrfA", "nrfH"),
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.nrf)

color.vector.nar <- c(cb.translator["vermillion"], cb.translator["yellow"], cb.translator["orange"], cb.translator["blue"])
names(color.vector.nar) <- c("narG", "narH", "narI", "napA")
nar.plot<- plot.scaffold.coverage(genesOfInterest = c("narG", "narH", "narI", "napA"),
                                  show.mean.coverage = FALSE,
                                  color.vector.to.use = color.vector.nar)

color.vector.nir <- c(cb.translator["vermillion"], cb.translator["yellow"])
names(color.vector.nir) <- c("nirS", "nirK")
nir.plot <- plot.scaffold.coverage(genesOfInterest = c("nirS", "nirK"),
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.nir)

color.vector.nox <- c(cb.translator["bluishgreen"])
names(color.vector.nir) <- c("nosZ")
nos.plot <- plot.scaffold.coverage(genesOfInterest = c("nosZ"),
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.nir)


pdf("results/2019_analysis_assembly/other_respiration/N_reduction_abundance.pdf",
    width = 10,
    height = 6)
(nrf.plot + nar.plot) / (nir.plot + nos.plot)
dev.off()
