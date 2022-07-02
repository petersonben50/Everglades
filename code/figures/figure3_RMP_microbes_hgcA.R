#### code/figures/figure3_RMP_microbes_hgcA.R ####
# Benjamin D. Peterson


#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ggpubr)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Data loading: RMP plot ####
RMP.hgcA <- readRDS("results/incubations/RMP_microbes_hgcA.rds")


#### Data loading: ordination plot ####
bc.ord.hgcA <- readRDS("results/metagenomes/assembly/hgcA/ordination.rds")


#### Data loading: metadata ####
metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
# Site metadata
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID
# Medium metadata
medium.metadata.vector <- metadata$medium
names(medium.metadata.vector) <- metadata$metagenomeID


#### Set up: color vector for microbes ####
color.code.df <- read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                           sheet = "color_codes")
color.code.vector.fused <- color.code.df$hexCode
names(color.code.vector.fused) <- color.code.df$clusterName
color.code.vector <- color.code.vector.fused[which(names(color.code.vector.fused) != "fused")]



#### Read in hgcA abundance and phylogenetic info ####
all.data <- read.csv("dataEdited/assembly_analysis/hgcA/depth/hgcA_coverage_scgNormalization.csv",
                           stringsAsFactors = FALSE) %>%
  full_join(read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                      sheet = "seq_class")) %>%
  arrange(classification) %>%
  filter(seqID %in% readLines("dataEdited/assembly_analysis/hgcA/hgcA_true.txt")) %>%
  filter(clusterName != "fused") %>%
  gather(key = MG,
         value = coverage,
         c(4:22)) %>%
  filter(MG != "total") %>%
  mutate(siteID = site.metadata.vector[MG],
         medium = medium.metadata.vector[MG]) %>%
  filter(medium == "sediment") %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  as.data.frame()





#### Plot: overall hgcA abundance ####
overall.hgcA.plot <- all.data %>%
  group_by(MG, siteID) %>%
  summarize(coverage = sum(coverage)) %>%
  group_by(siteID) %>%
  summarise(coverage_mean = mean(coverage),
            coverage_sd = sd(coverage),
            coverage_count = n(),
            coverage_se = coverage_sd / sqrt(coverage_count)) %>%
  ggplot(aes(y = coverage_mean,
             x = siteID)) + 
  geom_bar(fill = "grey75",
           stat="identity") +
  geom_errorbar(aes(ymin = coverage_mean - coverage_se,
                    ymax = coverage_mean + coverage_se),
                colour = "black",
                width = 0.33) +
  ylim(c(0, 16)) +
  ylab("hgcA gene abundance") +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none")



#### Plot: save out layout needed ####
figure.3.base.layer <- ggarrange(RMP.hgcA,
                                 ggarrange(bc.ord.hgcA + theme(legend.position = "none"),
                                           overall.hgcA.plot,
                                           ncol = 1,
                                           labels = c("B.", "C.")),
                                 ncol = 2,
                                 widths = c(1.5,1),
                                 labels = c("A.", ""))
pdf("results/figures/3/3_base_layer.pdf",
    width = 12,
    height = 7.5)
figure.3.base.layer
dev.off()


#### Data wrangling: Split "Unknown" category between "Unknown" and "Other" ####
unknown.seqs <- all.data %>%
  filter(grepl("Unknown", classification)) %>%
  select(seqID) %>%
  unlist(use.names = FALSE) %>%
  unique()
all.data[(all.data$clusterName =="Unknown") & !(all.data$seqID %in% unknown.seqs), "clusterName"] <- "Other"



#### Generate taxonomic bar plot ####
taxonomy.hgcA.plot <- all.data %>%
  mutate(clusterName = fct_relevel(clusterName,
                                   c("Aminicenantes", "Chloroflexi", "Firmicutes", "Spirochaetes",
                                     "Syntrophobacterales", "Methanogen",
                                     "Other", "Unknown"))) %>%
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
                    values = color.code.vector) + 
  ylim(c(0, 16)) +
  theme_classic() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = c(0.2, 0.8))
pdf("results/figures/3/3_taxonomy_layer.pdf",
    width = 4.8,
    height = 3.75)
taxonomy.hgcA.plot
dev.off()
