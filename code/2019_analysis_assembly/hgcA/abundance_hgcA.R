#### code/2019_analysis_assembly/hgcA/abundance_hgcA.R ####
# Benjamin D. Peterson

# Scripts to look at abundance of hgcA in a myriad of ways
# Can't believe I wrote "a myriad of ways"

#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(gridExtra)
library(lme4)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(tidyverse)



#### List of hgcA sequences ####

hgcA.list <- readLines("dataEdited/2019_analysis_assembly/hgcA/hgcA.txt")



#### Read in abundance and phylogenetic info ####

all.data <- read.csv("dataEdited/2019_analysis_assembly/hgcA/depth/hgcA_coverage.csv",
                               stringsAsFactors = FALSE) %>%
  full_join(read_xlsx("dataEdited/2019_analysis_assembly/hgcA/phylogeny/seq_classification.xlsx",
                      sheet = "seq_classification")) %>%
  arrange(classification) %>%
  filter(seqID %in% hgcA.list)






#### Read in median methylation fractions from incubations ####

inc.data <- readRDS("dataEdited/geochem/rel_methylation.rds")





#### Set up color vector ####

CB.color.vector <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")

color.code.df <- read_xlsx("dataEdited/2019_analysis_assembly/hgcA/phylogeny/seq_classification.xlsx",
                           sheet = "color_codes")
color.code.vector <- CB.color.vector[color.code.df$colorCode]
names(color.code.vector) <- color.code.df$clusterName



#### Read in metadata ####

metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")

# Site metadata
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID


# Medium metadata
medium.metadata.vector <- metadata$medium
names(medium.metadata.vector) <- metadata$metagenomeID



#### Generate vector of correct order of samples along sulfate gradient ####

MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")




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








#### Plot absolute abundance of each phylogenetic group at a site ####
# Plotting function for just overall abundance
all.data %>%
  group_by(MG, siteID, medium, classification) %>%
  summarize(coverage = round(sum(coverage), 4),
            count = n()) %>%
  ggplot(aes(x = MG,
             y = coverage,
             fill = classification,
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
    summarise(coverage = mean(coverage)) %>%
    ggplot(aes(y = coverage,
               x = siteID)) + 
    geom_bar(stat="identity") +
    labs(title = paste("hgcA abundance in ",
                       medium.of.interest,
                       sep = "")) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))
  }

hgcA.sediment <- plot.hgcA.overall.coverage(depth.data.of.interest = all.data,
                                            medium.of.interest = "sediment")
hgcA.porewater <- plot.hgcA.overall.coverage(depth.data.of.interest = all.data,
                                            medium.of.interest = "porewater")
hgcA.sediment / hgcA.porewater
pdf("results/2019_analysis_assembly/hgcA_abundance_overall.pdf",
      width = 5,
      height = 2.5)
hgcA.sediment
dev.off()



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
    geom_bar(stat="identity") +
    scale_fill_manual(name = "clusterName",
                      values = color.code.vector.of.interest) + 
    labs(title = paste("hgcA abundance in ",
                       medium.of.interest,
                       sep = "")) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))
  
}




hgcA.sediment <- plot.scaffold.coverage(depth.data.of.interest = all.data,
                                        medium.of.interest = "sediment")
hgcA.porewater <- plot.scaffold.coverage(depth.data.of.interest = all.data,
                                         medium.of.interest = "porewater")

pdf("results/2019_analysis_assembly/hgcA_abundance.pdf",
    width = 7.5,
    height = 4.5)
hgcA.sediment / hgcA.porewater
dev.off()









#### First focus on summary stats ####

# summary.stats.hgcA <- all.data %>%
#   group_by(siteID, medium, seqID) %>%
#   summarize(coverage = round(sum(coverage), 4),
#             count = n()) %>%
#   select(classification, count, total) %>%
#   mutate(`Percent of total\ncoverage` = round(total / sum(all.data$total) * 100, 2)) %>%
#   rename(`Taxonomic Group` = classification) %>%
#   rename(`Number of\nSequences` = count)

# # Write out a nice table with that data
# png("results/2019_analysis_assembly/hgcA/depth/taxonomic_cluster_summary_stats.png",
#     width = 5,
#     height = 2,
#     unit = "in",
#     res = 150)
# grid.table(summary.stats.hgcA,
#            rows = NULL)
# dev.off()











#### Plot abundance of hgcA in sediment against relative methylation ####

plotting.data.hgcA <- all.data %>%
  ungroup() %>%
  filter(medium == "sediment") %>%
  group_by(MG, siteID) %>%
  summarise(coverage = sum(coverage)) %>%
  group_by(siteID) %>%
  summarise(coverage_sd = sd(coverage),
            coverage = mean(coverage))

plotting.data.inc <- inc.data %>%
  group_by(siteID) %>%
  summarise(rel_meth_spike_sd = sd(rel_meth_spike),
            rel_meth_spike = median(rel_meth_spike))

plotting.data <- full_join(plotting.data.hgcA,
                           plotting.data.inc)


pdf("results/2019_analysis_assembly/relative_methylation_vs_hgcA_coverage_sediment.pdf",
    width = 7.5,
    height = 5)
plotting.data %>%
  ggplot(aes(x = coverage,
             y = rel_meth_spike,
             colour = siteID)) +
  geom_errorbar(aes(ymin = rel_meth_spike - rel_meth_spike_sd,
                    ymax = rel_meth_spike + rel_meth_spike_sd),
                colour = "black") +
  geom_errorbarh(aes(xmin = coverage - coverage_sd,
                     xmax = coverage + coverage_sd),
                colour = "black") +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "Read coverage of hgcA sequences\n(per 10 million reads)",
       y = "Median relative methylation percent\n(normalized to highest methylation %)",
       title = "Median relative methylation percent vs.\nhgcA coverage in sediment")
dev.off()





#### Run a linear mixed model ####


lmm.data <- all.data %>%
  filter(medium == "porewater") %>%
  group_by(MG, siteID) %>%
  summarise(coverage = sum(coverage)) %>%
  full_join(inc.data)

lmm.methylation <- lmer(meth_spike_per ~ coverage + (1|PW_source),
                        data = lmm.data)
summary(lmm.methylation)


plot(lmm.methylation)
qqnorm(resid(lmm.methylation))
qqline(resid(lmm.methylation)) 

