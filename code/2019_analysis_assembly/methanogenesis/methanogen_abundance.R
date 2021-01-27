#### code/2019_analysis_assembly/methanogenesis/methanogen_abundance.R ####
# Benjamin D. Peterson



#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(ggplot2)
library(patchwork)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Set site order for plotting ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")



#### Make plotting function ####

plot.scaffold.coverage <- function(scaffold.list,
                                   geneName,
                                   medium.of.interest) {
  depth.data %>%
    filter(scaffoldID %in% scaffold.list) %>%
    filter(medium == medium.of.interest) %>%
    mutate(siteID = fct_relevel(siteID, MG.order)) %>%
    group_by(siteID) %>%
    summarise(coverage = sum(coverage)) %>%
    ggplot(aes(y=coverage, x=siteID)) + 
    geom_bar(stat="identity") +
    labs(title = paste(geneName,
                       " in ",
                       medium.of.interest,
                       sep = "")) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))
  
}



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




#### Read in mcrA list ####

mcrA.scaffold.list <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/mcrA_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
# mtrA.scaffold.list <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/") %>%
#   strsplit("_[1-9]+") %>%
#   sapply("[", 1)



#### Generate coverage plots ####

mcrA.viz.sed <- plot.scaffold.coverage(mcrA.scaffold.list,
                                   medium.of.interest = "sediment",
                                   "mcrA")
mcrA.viz.pw <- plot.scaffold.coverage(mcrA.scaffold.list,
                                   medium.of.interest = "porewater",
                                   "mcrA")

mcrA.viz.sed / mcrA.viz.pw

#### Plot mcrA in sediment and porewater together ####

# png("results/2019_analysis_assembly/methanogenesis/mcrA_abundance.png",
#     units = 'in',
#     res = 200,
#     width = 5,
#     height = 5)
mcrA.viz.sed / mcrA.viz.pw
# dev.off()



#### Plot abundance with taxonomy ####

# Color vector first
color.vector.mcrA <- c(cb.translator["reddishpurple"],
                       cb.translator["orange"],
                       cb.translator["vermillion"])
names(color.vector.mcrA) <- read_xlsx("dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_clusters.xlsx") %>%
  select(clusterName) %>%
  unlist(use.names = FALSE) %>%
  unique()

depth.tax <- depth.data %>%
  filter(scaffoldID %in% mcrA.scaffold.list) %>%
  full_join(read_xlsx("dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_clusters.xlsx") %>%
              mutate(scaffoldID = seqID %>% strsplit("_[1-9]+") %>% sapply("[", 1))) %>%
  group_by(sampleID, siteID, medium, clusterName) %>%
  summarize(coverage = round(sum(coverage), 4)) %>%
  ungroup() %>%
  filter(medium == "sediment") %>%
  group_by(siteID, clusterName) %>%
  summarize(coverage = mean(coverage)) %>%
  ungroup() %>%
  mutate(siteID = fct_relevel(siteID, MG.order))

pdf("results/2019_analysis_assembly/methanogenesis/mcrA_abundance_tax.pdf",
    width = 6,
    height = 3)
depth.tax %>%
  ggplot(aes(x = siteID,
             y = coverage,
             fill = clusterName)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color.vector.mcrA) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Gene coverage\n(per 100X coverage of SCGs)",
       x = "Site ID",
       title = "Abundance of mcrA taxonomic clusters") +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black"))
dev.off()





#### Read in marker lists ####

list.of.marker.lists <- list.files(path = "dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/markers/",
                                   pattern = "methan_mark")
for (marker.of.interest in 1:length(list.of.marker.lists)) {
  
  scaffolds.of.interest = readLines(paste("dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/markers/",
                                          list.of.marker.lists[marker.of.interest],
                                          sep = "")) %>%
    strsplit("_[1-9]+") %>%
    sapply("[", 1)
  marker.df.temp <- data.frame(scaffoldID = scaffolds.of.interest,
                               geneName = rep(paste("marker_", marker.of.interest,
                                                    sep = ""),
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


#### Plot data ####
pdf("results/2019_analysis_assembly/methanogenesis/meth_markers_abundance.pdf",
    height = 4,
    width = 6)
marker.depth %>%
  filter(medium == "sediment") %>%
  group_by(scaffoldID, geneName, siteID) %>%
  summarise(coverage = mean(coverage)) %>%
  ungroup() %>%
  group_by(geneName, siteID) %>%
  summarise(coverage = sum(coverage)) %>%
  ungroup() %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  ggplot(aes(x = siteID,
             y = coverage)) +
  geom_line(mapping = aes(group = geneName)) +
  geom_point() +
  theme_classic() +
  stat_summary(geom = "point", fun = "mean",
               col = "black", fill = "red",
               size = 3, shape = 24) +
  labs(title = "Coverage of methanogenic marker genes",
       y = "Normalized to SCG coverage (100X)",
       x = "Site ID")
dev.off()



