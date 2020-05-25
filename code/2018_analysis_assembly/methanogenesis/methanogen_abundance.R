#### code/2018_analysis_assembly/methanogenesis/methanogen_abundance.R ####
# Benjamin D. Peterson




#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(ggplot2)
library(patchwork)
library(readxl)
library(tidyverse)






#### Set site order for plotting ####
MG.order <- c("WCA-2A-P", "WCA-2A-A", "WCA-3A-O", "WCA-3A-F", "WW")






#### Make plotting function ####

plot.scaffold.coverage <- function(scaffold.list,
                                   geneName) {

  depth.data %>%
    filter(scaffoldID %in% scaffold.list) %>%
    gather(key = sampleID,
           value = coverage,
           -1) %>%
    mutate(sampleID = fct_relevel(sampleID, MG.order)) %>%
    ggplot(aes(y=coverage, x=sampleID)) +
    geom_bar(position="stack", stat="identity") +
    labs(title = geneName) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))

}







##### Read in depth info ####

depth.data <- read.csv("dataEdited/2018_analysis_assembly/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE)





##### Read in MG metadata ####

MG.metadata <- read_xlsx("metadata/metagenomes/2018_MGs.xlsx")
renaming.vector <- MG.metadata$siteID
names(renaming.vector) <- MG.metadata$metagenomeID






##### Rename MGs with sample info ####

colnames(depth.data)[-1] <- renaming.vector[colnames(depth.data)[-1]]




#### Read in lists ####

mcrA.scaffold.list <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA/mcrA_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
mtrA.scaffold.list <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mtrA/mtrA_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)



#### Generate coverage plots ####

mcrA.viz <- plot.scaffold.coverage(mcrA.scaffold.list,
                                   "mcrA")
mtrA.viz <- plot.scaffold.coverage(mtrA.scaffold.list,
                                   "mtrA")

mcrA.viz / mtrA.viz

#### Plot mcrA and mtrA together ####

png("results/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA_mtrA_abundance.png",
    units = 'in',
    res = 200,
    width = 7.5,
    height = 3)
mcrA.viz / mtrA.viz
dev.off()








#### Read in marker lists ####

marker.list <- list()
plot.list <- list()
list.of.marker.lists <- list.files(path = "dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/markers/",
                                   pattern = "methan_mark")
for (marker.of.interest in 1:length(list.of.marker.lists)) {

         scaffolds.of.interest = readLines(paste("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/markers/",
                                                 list.of.marker.lists[marker.of.interest],
                                                 sep = "")) %>%
           strsplit("_[1-9]+") %>%
           sapply("[", 1)

         plot.list[[marker.of.interest]] <- plot.scaffold.coverage(scaffolds.of.interest,
                                                                   paste("Methanogenic marker ",
                                                                         marker.of.interest,
                                                                         sep = ""))
       }


png("results/2018_analysis_assembly/metabolicProteins/methanogenesis/meth_markers_abundance.png",
    units = "in",
    res = 200,
    height = 15,
    width = 12)
(plot.list[[1]] + plot.list[[2]]) /  (plot.list[[3]] + plot.list[[4]]) /  (plot.list[[5]] + plot.list[[6]]) /  (plot.list[[7]] + plot.list[[8]])
dev.off()
