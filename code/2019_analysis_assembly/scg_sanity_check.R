#### code/2019_analysis_assembly/scg_sanity_check.R ####
# Benjamin D. Peterson




#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(ggplot2)
library(patchwork)
library(readxl)
library(tidyverse)






#### Set site order for plotting ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")






#### Make plotting function ####

plot.scaffold.coverage <- function(scaffold.list,
                                   geneName,
                                   medium.of.interest) {
  testing <- depth.data %>%
    filter(scaffoldID %in% scaffold.list) %>%
    filter(medium == medium.of.interest) %>%
    mutate(siteID = fct_relevel(siteID, MG.order)) %>%
    group_by(sampleID, siteID) %>%
    summarise(coverage = sum(coverage)) %>%
    ungroup() %>%
    group_by(siteID) %>%
    summarise(coverage = mean(coverage)) %>%
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







#### Housekeeper gene plots ####

recA.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/houseKeeping/recA_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
recA.tree.viz <- plot.scaffold.coverage(recA.scaffolds,
                                        "recA",
                                        medium.of.interest = "sediment")
rpL14a.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/houseKeeping/rpL14a_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
rpL14b.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/houseKeeping/rpL14b_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
rpL14.tree.viz <- plot.scaffold.coverage(c(rpL14a.scaffolds, rpL14b.scaffolds),
                                         "rpL14",
                                         medium.of.interest = "sediment")

rpS3a.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/houseKeeping/rpS3a_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
rpS3b.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/houseKeeping/rpS3b_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
rpS3.tree.viz <- plot.scaffold.coverage(scaffold.list = c(rpS3a.scaffolds, rpS3b.scaffolds),
                                        geneName = "rpS3",
                                        medium.of.interest = "sediment")
pdf("results/2019_analysis_assembly/scg_sanity_check.pdf",
    width = 6,
    height = 9)
recA.tree.viz / rpL14.tree.viz / rpS3.tree.viz
dev.off()
