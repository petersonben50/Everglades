#### code/assembly_analysis/hgcA/hgcA_tree.R ####
# Benjamin D. Peterson

# This script will generate hgcA trees from the assembly-based
# hgcA analyses.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ggtree)
library(readxl)
library(tidyverse)
library(treeio)


#### Read in tree file ####
hgcA.tree <- readRDS("dataEdited/assembly_analysis/hgcA/phylogeny/final/hgcA_clean_tree.rds")


#### Read in color vector ####
color.vector <- readRDS("dataEdited/assembly_analysis/hgcA/phylogeny/final/hgcA_clean_tree_color_vector.rds")


#### List of hgcA sequences ####
hgcA.list <- readLines("dataEdited/assembly_analysis/hgcA/hgcA.txt")


#### Read in binning information ####
hgcA2binTax.df <- read.table("/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/bin_analysis/binning_initial/binsFinal/hgcA_bin_key.tsv",
                          sep = '\t',
                          header = TRUE) %>%
  select(binID, hgcA) %>%
  left_join(readRDS("dataEdited/assembly_analysis/hgcA/dereplication/hgcA_derep_key.rds")) %>%
  left_join(read_xlsx("dataEdited/bin_analysis/binning_initial/binsGood/taxonomy_summary.xlsx")) %>%
  select(hgcA_rep, taxID)
hgcA2binTax.vector <- paste("Binned: ",
                            hgcA2binTax.df$taxID,
                            sep = "")
names(hgcA2binTax.vector) <- hgcA2binTax.df$hgcA_rep


#### Replace names of binned seqs ####
binned.index <- which(hgcA.tree$tip.label %in% names(hgcA2binTax.vector))
names(hgcA2binTax.vector)[which(!(names(hgcA2binTax.vector) %in% hgcA.tree$tip.label[binned.index]))]
hgcA.tree$tip.label[binned.index] <- paste(hgcA.tree$tip.label[binned.index],
                                           " (",
                                           hgcA2binTax.vector[hgcA.tree$tip.label[binned.index]],
                                           ")",
                                           sep = "")


#### Make tree ####
hgcA.tree.plot <- ggtree(hgcA.tree,
                         aes(x = 0,
                             xend = 5)) + 
  geom_tiplab(size=1.5, align = TRUE, colour = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2) +
  geom_treescale(x = 0.05,
                 y = 120,
                 width = 0.5)
pdf("results/metagenomes/assembly/hgcA/hgcA_tree_RAxML_rooted_with_bins.pdf",
    height = 60,
    width = 10)
hgcA.tree.plot
dev.off()



#### Read in and prep depth info ####
depth.data <- read.csv("dataEdited/assembly_analysis/hgcA/depth/hgcA_coverage_scgNormalization.csv",
                       stringsAsFactors = FALSE)


#### Add node labels to depth information ####
hgcA.node.labels <- grep(paste(hgcA.list,
                               collapse="|"),
                         hgcA.tree$tip.label,
                         value=TRUE)
names(hgcA.node.labels) <- hgcA.node.labels %>%
  strsplit(" ") %>% sapply("[", 1) 
names(hgcA.node.labels) <- paste(strsplit(names(hgcA.node.labels),
                                          split = "_",) %>% 
                                   sapply("[", 1),
                                 strsplit(names(hgcA.node.labels),
                                          split = "_",) %>% 
                                   sapply("[", 2),
                                 sep = "_")

depth.data.hgcA <- depth.data %>%
  filter(scaffoldID %in% names(hgcA.node.labels))
rownames(depth.data.hgcA) <- hgcA.node.labels[depth.data.hgcA$scaffoldID]
depth.data.hgcA <- depth.data.hgcA %>% select(-c(scaffoldID, seqID, length, total))


#### Rename MG with site ID ####
MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
renaming.vector <- paste(MG.metadata$siteID,
                         MG.metadata$metagenomeID,
                         sep = "_")
names(renaming.vector) <- MG.metadata$metagenomeID
colnames(depth.data.hgcA) <- renaming.vector[colnames(depth.data.hgcA)]
depth.data.hgcA <- depth.data.hgcA %>%
  select(`2A-N_KMBP005A`, `2A-N_KMBP005G`,
         `2A-A_KMBP005B`, `2A-A_KMBP005H`,
         `3A-O_KMBP005F`, `3A-O_KMBP005L`,
         `3A-N_KMBP005E`, `3A-N_KMBP005K`,
         `3A-F_KMBP005D`, `3A-F_KMBP005J`,
         `LOX8_KMBP005C`, `LOX8_KMBP005I`)


#### Generate tree with depth ####
pdf("results/metagenomes/assembly/hgcA/hgcA_tree_RAxML_rooted_with_bins_depth.pdf",
    height = 10,
    width = 7.5)
gheatmap(hgcA.tree.plot,
         depth.data.hgcA, offset=1, width=0.5, font.size=2, 
         colnames_angle=-90, hjust=0)  +
  scale_fill_viridis_c(option="D", name="continuous\nvalue")
dev.off()


#### Generate tree with relative depth ####
rel.depth.data.hgcA <- depth.data.hgcA / rowSums(depth.data.hgcA)
pdf("results/metagenomes/assembly/hgcA/hgcA_tree_RAxML_rooted_with_bins_depthRelative.pdf",
    height = 10,
    width = 7.5)
gheatmap(hgcA.tree.plot,
         rel.depth.data.hgcA, offset=1, width=0.5, font.size=2, 
         colnames_angle=-90, hjust=0)  +
  scale_fill_viridis_c(option="D", name="continuous\nvalue")
dev.off()

