#### code/2018_analysis_assembly/sulfur/SRB_abundance.R ####
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





#### Read in dsr lists ####

dsrA.scaffolds <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_red_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)

dsrD.scaffolds <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrD_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)



#### Plot out dsr depths ####

dsrA.tree.viz <- plot.scaffold.coverage(dsrA.scaffolds,
                                        "dsrA")
dsrD.tree.viz <- plot.scaffold.coverage(dsrD.scaffolds,
                                        "dsrD")




#### Plot dsrA and dsrD together ####

pdf("results/2018_analysis_assembly/metabolicProteins/sulfur/dsr_reductive_abund.pdf",
    width = 7.5,
    height = 3)
dsrA.tree.viz + dsrD.tree.viz
dev.off()





#### Read in lists of asr subunits ####

asrA.scaffolds <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/asr/sulfite_red_A_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
asrB.scaffolds <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/asr/sulfite_red_B_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
asrC.scaffolds <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/asr/sulfite_red_C_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)



