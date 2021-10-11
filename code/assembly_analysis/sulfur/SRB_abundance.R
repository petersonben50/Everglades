#### code/assembly_analysis/sulfur/SRB_abundance.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(ggplot2)
library(patchwork)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")
source("code/assembly_analysis/metabolic_protein_plots.R")


##### Read in MG metadata ####
MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.renaming.vector <- MG.metadata$siteID
names(site.renaming.vector) <- MG.metadata$metagenomeID
medium.renaming.vector <- MG.metadata$medium
names(medium.renaming.vector) <- MG.metadata$metagenomeID


##### Read in depth info ####
depth.data <- read.csv("dataEdited/assembly_analysis/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE) %>%
  gather(key = sampleID,
         value = coverage,
         -1) %>%
  mutate(siteID = site.renaming.vector[sampleID]) %>%
  mutate(medium = medium.renaming.vector[sampleID])


#### Read in dsr lists ####
dsrA.scaffolds <- readLines("dataEdited/assembly_analysis/metabolicProteins/sulfur/dsr/dsrA/dsrA_red_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)
dsrD.scaffolds <- readLines("dataEdited/assembly_analysis/metabolicProteins/sulfur/dsr/dsrD_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)


#### Read in gene lists ####
list.of.marker.lists <- c("dataEdited/assembly_analysis/metabolicProteins/sulfur/dsr/dsrA/dsrA_red_list.txt",
                          "dataEdited/assembly_analysis/metabolicProteins/sulfur/dsr/dsrD_derep_list.txt")
names(list.of.marker.lists) <- c("dsrA", "dsrD")

for (marker.of.interest in 1:length(list.of.marker.lists)) {
  
  scaffolds.of.interest = readLines(list.of.marker.lists[marker.of.interest]) %>%
    strsplit("_[1-9]+") %>%
    sapply("[", 1)
  gene_name <- names(list.of.marker.lists)[marker.of.interest]
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



#### Plot out dsr depths ####
# dsrA.tree.viz <- plot.scaffold.coverage(dsrA.scaffolds,
#                                         "dsrA",
#                                         medium.of.interest = "sediment")
# dsrD.tree.viz <- plot.scaffold.coverage(scaffold.list = dsrD.scaffolds,
#                                         geneName = "dsrD",
#                                         medium.of.interest = "sediment")


#### dsrD plots ####
color.vector.dsrD <- c(cb.translator["blue"])
names(color.vector.dsrD) <- c("dsrD")
dsrD.plot <- plot.scaffold.coverage(genesOfInterest = names(color.vector.dsrD),
                                   show.mean.coverage = FALSE,
                                   color.vector.to.use = color.vector.dsrD)


#### Plot dsrA and dsrD together ####
# pdf("results/assembly_analysis/sulfur/dsrD_abund.pdf",
#     width = 5,
#     height = 3)
dsrD.plot
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


#### Plot dsrA and dsrD together ####
dsr.color.vector <- cb.translator[c("blue", "skyblue")]
names(dsr.color.vector) <- c("dsrA", "dsrD")
dsr.plots <- dsr.abundance %>%
  ggplot(aes(y = coverage,
             x = siteID,
             fill = geneName)) +
  geom_bar(stat="identity",
           position = "dodge") +
  scale_fill_manual(values = dsr.color.vector,
                    name = "Gene ID") +
  theme_bw() +
  labs(y = "Gene coverage\n(per 100X coverage of SCGs)") +
  theme(axis.text.x = element_text(colour="black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour="black"),
        legend.position = c(0.5, 0.7))

# Save out plot
saveRDS(dsr.plots,
        "results/metagenomes/assembly/SRB/dsr_abundance.rds")




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
# 
# # Color vector first
# color.vector.dsrA <- c(cb.translator["skyblue"],
#                        cb.translator["blue"],
#                        cb.translator["bluishgreen"],
#                        cb.translator["black"])
# names(color.vector.dsrA) <- c("cluster1", "cluster2", "cluster3", "cluster4")
# 
# depth.tax <- depth.data %>%
#   filter(scaffoldID %in% dsrA.scaffolds) %>%
#   filter(medium == "sediment") %>%
#   full_join(read.csv("dataEdited/assembly_analysis/metabolicProteins/sulfur/dsr/reductive_dsrA_phylogeny/dsrA_cluster.csv",
#                      stringsAsFactors = FALSE) %>%
#               mutate(scaffoldID = seqID %>% strsplit("_[1-9]+") %>% sapply("[", 1))) %>%
#   group_by(sampleID, siteID, clusterID) %>%
#   summarize(coverage = round(sum(coverage), 4)) %>%
#   ungroup() %>%
#   group_by(siteID, clusterID) %>%
#   summarize(coverage = mean(coverage)) %>%
#   ungroup() %>%
#   mutate(siteID = fct_relevel(siteID, MG.order))
# 
# pdf("results/assembly_analysis/sulfur/dsrA_abundance_tax.pdf",
#     width = 6,
#     height = 3)
# depth.tax %>%
#   ggplot(aes(x = siteID,
#              y = coverage,
#              fill = clusterID)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = color.vector.dsrA) +
#   theme_classic() +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(y = "Gene coverage\n(per 100X coverage of SCGs)",
#        x = "Site ID",
#        title = "Abundance of dsrA phylogenetic clusters") +
#   theme(axis.text.x = element_text(colour = "black")) +
#   theme(axis.text.y = element_text(colour = "black"))
# dev.off()
# 
