#### code/2018_analysis_assembly/hgcA/hgcA_phylogeny_abundance.R ####
# Benjamin D. Peterson

# This script will generate a phylogenetic tree with the 
# corresponding abundance information in it. 
##########################


#### Always start with a clean slate ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(ggtree)
library(phangorn)
library(tidyverse)
library(treeio)



#### Read in needed data ####

# Read in tree:
hgcA.tree <- readRDS("dataEdited/2018_analysis_assembly/hgcA/phylogeny/hgcA_clean_tree.rds")

# Read in color vector
color.vector <- readRDS("dataEdited/2018_analysis_assembly/hgcA/phylogeny/hgcA_clean_tree_color_vector.rds")

# Read in depth data
depth.data <- read.csv("dataEdited/2018_analysis_assembly/hgcA/depth/hgcA_coverage_final.csv",
                       stringsAsFactors = FALSE) %>%
  select(-length)

# Read in list of final hgcA seqs
hgcA.list <- readLines("dataEdited/2018_analysis_assembly/hgcA/hgcA.txt")






#### Clean up depth data ####

# Make into dataframe for plotting
depth.for.plotting <- data.frame(id = depth.data$seqID,
                                 value = depth.data$total,
                                 stringsAsFactors = FALSE)




#### Visualize tree ####
# Visualize just the tree
tree <- ggtree(hgcA.tree, aes(x = 0,
                              xend = 3.5)) + 
  geom_tiplab(size=2, align = TRUE, col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
# Add in depth
treeAndDepth <- facet_plot(tree,
                           panel = "Coverage",
                           data = depth.for.plotting,
                           geom=geom_segment, aes(x = 0,
                                                  xend = value,
                                                  y = y,
                                                  yend = y),
                           size = 1) + theme_tree2()

png("results/2018_analysis_assembly/hgcA/phylogeny/hgcA_tree_abundance.png",
    units = "in",
    res = 240,
    height = 10,
    width = 6)
treeAndDepth
dev.off()
