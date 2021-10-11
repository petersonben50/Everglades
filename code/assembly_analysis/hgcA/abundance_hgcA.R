#### code/assembly_analysis/hgcA/abundance_hgcA.R ####
# Benjamin D. Peterson

# Scripts to look at abundance of hgcA in a myriad of ways
# Can't believe I wrote "a myriad of ways"

#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### List of hgcA sequences ####
hgcA.list <- readLines("dataEdited/assembly_analysis/hgcA/hgcA_true.txt")


#### Read in abundance and phylogenetic info ####
all.data.fused <- read.csv("dataEdited/assembly_analysis/hgcA/depth/hgcA_coverage_scgNormalization.csv",
                               stringsAsFactors = FALSE) %>%
  full_join(read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                      sheet = "seq_class")) %>%
  arrange(classification) %>%
  filter(seqID %in% hgcA.list) 


#### Set up color vector ####
CB.color.vector <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
color.code.df <- read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                           sheet = "color_codes")
color.code.vector.fused <- CB.color.vector[color.code.df$colorCode]
names(color.code.vector.fused) <- color.code.df$clusterName



#### Read in metadata ####
metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
# Site metadata
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID
# Medium metadata
medium.metadata.vector <- metadata$medium
names(medium.metadata.vector) <- metadata$metagenomeID



#### Check on abundance fused hgcAB ####
all.data.fused %>%
  filter(clusterName == "fused") %>%
  gather(key = MG,
         value = coverage,
         c(4:22)) %>%
  group_by(MG) %>%
  summarise(coverage = sum(coverage)) %>%
  mutate(siteID = site.metadata.vector[MG],
         medium = medium.metadata.vector[MG]) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  filter(medium == "sediment") %>%
  ggplot(aes(x = siteID,
             y = coverage)) +
  geom_point() +
  theme_classic()


#### Remove fused sequences ####
all.data <- all.data.fused %>%
  filter(clusterName != "fused")
color.code.vector <- color.code.vector.fused[unique(all.data$clusterName)]


#### Prepare data for plotting ####
all.data <- all.data %>%
  gather(key = MG,
         value = coverage,
         c(4:22)) %>%
  filter(MG != "total") %>%
  mutate(siteID = site.metadata.vector[MG],
         medium = medium.metadata.vector[MG]) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  as.data.frame()


#### Save out R object of average coverage across sites for use with incubations ####
all.data %>%
  filter(medium == "sediment") %>%
  group_by(MG, siteID) %>%
  summarise(coverage = sum(coverage)) %>%
  ungroup() %>%
  saveRDS("dataEdited/assembly_analysis/hgcA/hgcA_abundance_site.rds")



#### Plot absolute abundance of each phylogenetic group at a site ####
# Plotting function for just overall abundance
all.data %>%
  group_by(MG, siteID, medium, clusterName) %>%
  summarize(coverage = round(sum(coverage), 4),
            count = n()) %>%
  ggplot(aes(x = MG,
             y = coverage,
             fill = clusterName,
             width = 3)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~siteID, nrow = 1) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# Replicates look fairly similar, so let's make bar plots
# where we average the two replicates together.



#### Plot out total hgcA abundance ####

plot.hgcA.overall.coverage <- function(depth.data.of.interest = all.data,
                                       medium.of.interest = "sediment") {
  depth.data.of.interest %>%
    mutate(clusterName = fct_relevel(clusterName,
                                     names(color.code.vector)[length(names(color.code.vector)):1])) %>%
    filter(medium == medium.of.interest) %>%
    group_by(MG, siteID) %>%
    summarize(coverage = sum(coverage)) %>%
    group_by(siteID) %>%
    summarise(coverage_mean = mean(coverage),
              coverage_sd = sd(coverage),
              coverage_count = n(),
              coverage_se = coverage_sd / sqrt(coverage_count)) %>%
    ggplot(aes(y = coverage_mean,
               x = siteID)) + 
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin = coverage_mean - coverage_se,
                      ymax = coverage_mean + coverage_se),
                  colour = "black",
                  width = 0.33) +
    labs(title = paste("hgcA abundance in ",
                       medium.of.interest,
                       sep = "")) +
    ylim(c(0, 16)) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"),
          legend.position = "none")
  }

hgcA.sediment <- plot.hgcA.overall.coverage(depth.data.of.interest = all.data,
                                            medium.of.interest = "sediment")
hgcA.porewater <- plot.hgcA.overall.coverage(depth.data.of.interest = all.data,
                                            medium.of.interest = "porewater")
hgcA.sediment / hgcA.porewater



#### Plotting function with taxa ####
plot.scaffold.coverage <- function(depth.data.of.interest = all.data,
                                   medium.of.interest = "sediment",
                                   color.code.vector.of.interest = color.code.vector) {
  depth.data.of.interest %>%
    mutate(clusterName = fct_relevel(clusterName,
                                     names(color.code.vector)[length(names(color.code.vector)):1])) %>%
    filter(medium == medium.of.interest) %>%
    group_by(MG, siteID, clusterName) %>%
    summarize(coverage = sum(coverage)) %>%
    group_by(siteID, clusterName) %>%
    summarise(coverage = mean(coverage)) %>%
    ggplot(aes(y = coverage,
               x = siteID,
               fill = clusterName)) + 
    geom_bar(stat = "identity",
             position = "dodge") +
    scale_fill_manual(name = "clusterName",
                      values = color.code.vector.of.interest) + 
    labs(title = paste("hgcA abundance in ",
                       medium.of.interest,
                       sep = "")) +
    ylim(c(0, 7)) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"),
          legend.position = "none")
  
}




hgcA.sediment <- plot.scaffold.coverage(depth.data.of.interest = all.data,
                                        medium.of.interest = "sediment")
hgcA.porewater <- plot.scaffold.coverage(depth.data.of.interest = all.data,
                                         medium.of.interest = "porewater")

# pdf("results/metagenomes/assembly/hgcA/hgcA_abundance_sediment.pdf",
#     width = 5,
#     height = 3.5)
hgcA.sediment
# dev.off()



#### Relative coverage of each group ####

tax.group.coverage <- all.data %>%
  mutate(clusterName = fct_relevel(clusterName,
                                   names(color.code.vector)[length(names(color.code.vector)):1])) %>%
  filter(medium == "sediment") %>%
  group_by(MG, siteID, clusterName) %>%
  summarize(coverage = sum(coverage)) %>%
  group_by(siteID, clusterName) %>%
  summarise(coverage = mean(coverage)) %>%
  spread(key = siteID,
         value = coverage) %>%
  mutate(total = `2A-N` + `2A-A` + `3A-O` +
           `3A-N` + `3A-F` + `LOX8`) %>%
  arrange(desc(total))

tax.group.coverage.relative <- tax.group.coverage
tax.group.coverage.relative[, -1] <- sapply(names(tax.group.coverage.relative)[-1],
                                            function(siteID) {
                                              unlist(tax.group.coverage.relative[, siteID] / colSums(tax.group.coverage.relative[, -1])[siteID],
                                                     use.names = FALSE)
                                            })



#### Plot relative abundance of each hgcA group over the gradient ####
tax.group.coverage.relative.plot <- tax.group.coverage.relative %>%
  gather(key = siteID,
         value = rel.coverage,
         -1) %>%
  filter(siteID != "total") %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  ggplot(aes (x = siteID,
              y = rel.coverage,
              group = clusterName)) +
  geom_point(aes(color = clusterName)) + 
  geom_line(aes(color = clusterName)) +
  scale_color_manual(name = "hgcA classification",
                    values = color.code.vector) +
  theme_classic() +
  ylim(0, 0.6) +
  labs(y = "hgcA relative abundance") +
  theme(axis.text.x = element_text(colour="black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour="black"))
saveRDS(tax.group.coverage.relative.plot,
        "results/metagenomes/assembly/hgcA/tax_group_rel_coverage.rds")
