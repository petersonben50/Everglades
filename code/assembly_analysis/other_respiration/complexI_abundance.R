#### code/2019_analysis_assembly/other/complexI_abundance.R ####
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


# 
# #### Make plotting function ####
# 
# plot.scaffold.coverage <- function(scaffold.list,
#                                    geneName,
#                                    medium.of.interest,
#                                    use.points = FALSE) {
#   temp.data <- depth.data %>%
#     filter(scaffoldID %in% scaffold.list) %>%
#     filter(medium == medium.of.interest) %>%
#     mutate(siteID = fct_relevel(siteID, MG.order))
#   if (use.points == TRUE) {
#     temp.data %>%
#       group_by(sampleID, siteID) %>%
#       summarise(coverage = sum(coverage)) %>%
#       ggplot(aes(y=coverage, x=siteID)) + 
#       geom_point(stat="identity") +
#       labs(title = paste(geneName,
#                          " in ",
#                          medium.of.interest,
#                          sep = "")) +
#       theme_classic() +
#       theme(axis.text.x = element_text(colour="black"),
#             axis.text.y = element_text(colour="black"))
#   } else {
#     temp.data %>%
#       group_by(siteID) %>%
#       summarise(coverage = sum(coverage)) %>%
#       ggplot(aes(y=coverage, x=siteID)) + 
#       geom_bar(stat="identity") +
#       labs(title = paste(geneName,
#                          " in ",
#                          medium.of.interest,
#                          sep = "")) +
#       theme_classic() +
#       theme(axis.text.x = element_text(colour="black"),
#             axis.text.y = element_text(colour="black"))
#   }
# }
# 
# 






##### Read in MG metadata ####

MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.renaming.vector <- MG.metadata$siteID
names(site.renaming.vector) <- MG.metadata$metagenomeID

medium.renaming.vector <- MG.metadata$medium
names(medium.renaming.vector) <- MG.metadata$metagenomeID





#### Read in gene lists ####
list.of.marker.lists <- list.files(path = "dataEdited/2019_analysis_assembly/metabolicProteins/other_respiration/complexI/",
                                   pattern = "_derep_list.txt")
for (marker.of.interest in 1:length(list.of.marker.lists)) {
  
  scaffolds.of.interest = readLines(paste("dataEdited/2019_analysis_assembly/metabolicProteins/other_respiration/complexI/",
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




##### Read in depth info ####

depth.data <- read.csv("dataEdited/2019_analysis_assembly/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE) %>%
  gather(key = sampleID,
         value = coverage,
         -1) %>%
  mutate(siteID = site.renaming.vector[sampleID]) %>%
  mutate(medium = medium.renaming.vector[sampleID])






#### nuo plots ####
color.vector.nuo <- c(cb.translator["bluishgreen"], cb.translator["blue"], cb.translator["skyblue"], cb.translator["black"])
names(color.vector.nuo) <- c("nuoB", "nuoC", "nuoD", "nuoI")
nuo.B.plot <- plot.scaffold.coverage(genesOfInterest =  c("nuoB", "nuoC", "nuoD", "nuoI"),
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.nuo)

color.vector.nuo <- c(cb.translator["bluishgreen"], cb.translator["blue"], cb.translator["skyblue"])
names(color.vector.nuo) <- c("nuoE", "nuoF", "nuoG")
nuo.plot <- plot.scaffold.coverage(genesOfInterest = c("nuoE", "nuoF", "nuoG"),
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.nuo)

color.vector.nuo <- c(cb.translator["bluishgreen"], cb.translator["blue"], cb.translator["skyblue"])
names(color.vector.nuo) <- c("nuoL", "nuoM", "nuoN")
nuo.M.plot <- plot.scaffold.coverage(genesOfInterest = c("nuoL", "nuoM", "nuoN"),
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.nuo)

pdf("results/2019_analysis_assembly/other_respiration/complexI.pdf",
    width = 6,
    height = 8)
nuo.B.plot / nuo.plot / nuo.M.plot
dev.off()

