#### code/2019_analysis_assembly/sulfur/SRB_abundance.R ####
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


#### Read in dsr lists ####
dsrA.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/dsrA_red_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
dsrD.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrD_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)



#### Plot out dsr depths ####
dsrA.tree.viz <- plot.scaffold.coverage(dsrA.scaffolds,
                                        "dsrA",
                                        medium.of.interest = "sediment")
dsrD.tree.viz <- plot.scaffold.coverage(scaffold.list = dsrD.scaffolds,
                                        geneName = "dsrD",
                                        medium.of.interest = "sediment")


#### Plot dsrA and dsrD together ####
# pdf("results/2019_analysis_assembly/sulfur/dsr_reductive_abund.pdf",
#     width = 5,
#     height = 5)
dsrA.tree.viz / dsrD.tree.viz
# dev.off()


#### Plot dsrA and dsrD side-by-side ####
dsrA.df <- data.frame(scaffoldID = dsrA.scaffolds,
                      geneName = rep("dsrA",
                                     length(dsrA.scaffolds)))
dsrD.df <- data.frame(scaffoldID = dsrD.scaffolds,
                      geneName = rep("dsrD",
                                     length(dsrD.scaffolds)))
dsr.df <- rbind(dsrA.df, dsrD.df)
dsr.abundance <- left_join(dsr.df,
                           depth.data) %>%
  mutate(siteID = site.renaming.vector[sampleID]) %>%
  mutate(medium = medium.renaming.vector[sampleID]) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  filter(medium == "sediment") %>%
  group_by(sampleID, siteID, medium, geneName) %>%
  summarise(coverage = sum(coverage)) %>%
  ungroup() %>%
  group_by(siteID, geneName) %>%
  summarise(coverage = mean(coverage))


# pdf("results/2019_analysis_assembly/sulfur/dsr_reductive_abund_sideBySide.pdf",
#     width = 5,
#     height = 3)
dsr.abundance %>%
  ggplot(aes(y = coverage, x = siteID,
             fill = geneName)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_bw()
# dev.off()  



#### Compare dsrA to sulfide ####
read_xlsx("dataEdited/geochem/geochem_data_2019.xlsx",
                             sheet = "PW_SW_geochem") %>%
  filter(medium == "PW") %>%
  select(siteID, Sulfide_µg.L, Sulfate_mg.L) %>%
  left_join(dsr.abundance) %>%
  filter(geneName == "dsrA") %>%
  ggplot(aes(x = coverage,
             y = log(Sulfide_µg.L, 10))) +
  geom_point() +
  theme_bw()



#### Plot dsrA with phylogenetic info ####

# Color vector first
color.vector.dsrA <- c(cb.translator["skyblue"],
                       cb.translator["blue"],
                       cb.translator["bluishgreen"],
                       cb.translator["black"])
names(color.vector.dsrA) <- c("cluster1", "cluster2", "cluster3", "cluster4")

depth.tax <- depth.data %>%
  filter(scaffoldID %in% dsrA.scaffolds) %>%
  filter(medium == "sediment") %>%
  full_join(read.csv("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/reductive_dsrA_phylogeny/dsrA_cluster.csv",
                     stringsAsFactors = FALSE) %>%
              mutate(scaffoldID = seqID %>% strsplit("_[1-9]+") %>% sapply("[", 1))) %>%
  group_by(sampleID, siteID, clusterID) %>%
  summarize(coverage = round(sum(coverage), 4)) %>%
  ungroup() %>%
  group_by(siteID, clusterID) %>%
  summarize(coverage = mean(coverage)) %>%
  ungroup() %>%
  mutate(siteID = fct_relevel(siteID, MG.order))

pdf("results/2019_analysis_assembly/sulfur/dsrA_abundance_tax.pdf",
    width = 6,
    height = 3)
depth.tax %>%
  ggplot(aes(x = siteID,
             y = coverage,
             fill = clusterID)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color.vector.dsrA) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Gene coverage\n(per 100X coverage of SCGs)",
       x = "Site ID",
       title = "Abundance of dsrA phylogenetic clusters") +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black"))
dev.off()



#### Read in lists of asr subunits ####
asrA.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/asr/sulfite_red_A_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
asrB.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/asr/sulfite_red_B_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
asrC.scaffolds <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/asr/sulfite_red_C_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)


asrA.tree.viz <- plot.scaffold.coverage(asrA.scaffolds,
                                        "asrA",
                                        medium.of.interest = "sediment")
asrB.tree.viz <- plot.scaffold.coverage(scaffold.list = asrB.scaffolds,
                                        geneName = "asrB",
                                        medium.of.interest = "sediment")
asrC.tree.viz <- plot.scaffold.coverage(scaffold.list = asrC.scaffolds,
                                        geneName = "asrC",
                                        medium.of.interest = "sediment")
asrA.tree.viz / asrB.tree.viz / asrC.tree.viz
