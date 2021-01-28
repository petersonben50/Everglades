#### code/2019_analysis_assembly/other/hydrogenase_abundance.R ####
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
depth.data <- read.csv("dataEdited/2019_analysis_assembly/metabolicProteins/hydrogenases/depth/hydrogenase_depth_clean.csv",
                       stringsAsFactors = FALSE) %>%
  gather(key = sampleID,
         value = coverage,
         -1) %>%
  mutate(siteID = site.renaming.vector[sampleID]) %>%
  mutate(medium = medium.renaming.vector[sampleID])





#### Read in marker lists ####

list.of.marker.lists <- list.files(path = "dataEdited/2019_analysis_assembly/metabolicProteins/hydrogenases/lists",
                                   pattern = "_derep_list.txt")
for (marker.of.interest in 1:length(list.of.marker.lists)) {
  scaffolds.of.interest = readLines(paste("dataEdited/2019_analysis_assembly/metabolicProteins/hydrogenases/lists/",
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



#### Hydrogenase classes ####
class.list <- read_xlsx("dataEdited/2019_analysis_assembly/metabolicProteins/hydrogenases/hydrogenase_classes.xlsx")

h2.evolvers <- class.list %>% filter(functionID == "evolution") %>% select(geneName) %>% unlist(use.names = FALSE) %>% gsub("_group_", "_", .)
h2.takerUppers <- class.list %>% filter(functionID == "uptake") %>% select(geneName) %>% unlist(use.names = FALSE) %>% gsub("_group_", "_", .)
h2.sensing <- class.list %>% filter(functionID == "sensing") %>% select(geneName) %>% unlist(use.names = FALSE) %>% gsub("_group_", "_", .)
h2.bi <- class.list %>% filter(functionID == "bidirectional") %>% select(geneName) %>% unlist(use.names = FALSE) %>% gsub("_group_", "_", .)

#### Generate plots ####

# H2-evolving
color.vector.H2E <- cb.translator[1:length(h2.evolvers)]
names(color.vector.H2E) <- h2.evolvers
H2E.plot <- plot.scaffold.coverage(marker.depth.df = marker.depth,
                                   genesOfInterest = h2.evolvers,
                                   medium.of.interest = "sediment",
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.H2E,
                                   MG.order.vector = MG.order,
                                   titleToUse = paste("H2-evolving hydrogenases"))

# H2-uptake                  
color.vector.H2U <- cb.translator[1:length(h2.takerUppers)]
names(color.vector.H2U) <- h2.takerUppers
H2U.plot <- plot.scaffold.coverage(marker.depth.df = marker.depth,
                                   genesOfInterest = h2.takerUppers,
                                   medium.of.interest = "sediment",
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.H2U,
                                   MG.order.vector = MG.order,
                                   titleToUse = paste("H2-consuming hydrogenases"))

# H2-sensing                  
color.vector.H2S <- cb.translator[1:length(h2.sensing)]
names(color.vector.H2S) <- h2.sensing
H2S.plot <- plot.scaffold.coverage(marker.depth.df = marker.depth,
                                   genesOfInterest = h2.sensing,
                                   medium.of.interest = "sediment",
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.H2S,
                                   MG.order.vector = MG.order,
                                   titleToUse = paste("H2-sensing hydrogenases"))

# H2-bidirectional                  
color.vector.H2B <- cb.translator[1:length(h2.bi)]
names(color.vector.H2B) <- h2.bi
H2B.plot <- plot.scaffold.coverage(marker.depth.df = marker.depth,
                                   genesOfInterest = h2.bi,
                                   medium.of.interest = "sediment",
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.H2B,
                                   MG.order.vector = MG.order,
                                   titleToUse = paste("H2-sensing hydrogenases"))


pdf("results/2019_analysis_assembly/non_respiration/hydrogenases_abundance.pdf",
    width = 10,
    height = 6)
(H2E.plot + H2S.plot) / (H2U.plot + H2B.plot)
dev.off()
