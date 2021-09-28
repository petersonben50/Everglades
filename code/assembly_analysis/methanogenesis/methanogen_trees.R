#### code/2019_analysis_assembly/methanogenesis/methanogen_trees.R ####
# Benjamin D. Peterson



#### Clean up on aisle R ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(ggtree)
library(phangorn)
library(tidyverse)
library(treeio)


#### Read in mcrA tree metadata ####



#### Read in mcrA tree ####

mcrA.tree.unrooted <- read.newick("dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_for_phylogeny.tree")
mcrA.tree <- midpoint(mcrA.tree.unrooted)


#### Read in mcrA list from this study ####

mcrA.list <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/mcrA_derep_list.txt")


#### Make indices ####

my.seq.index <- which(mcrA.tree$tip.label %in% mcrA.list)


#### Make color vector ####

color.vector <- rep("black", length(mcrA.tree$tip.label))
color.vector[my.seq.index] <- "red"



#### Make mcrA tree ####


mcrA.tree.object <- ggtree(mcrA.tree) +
  geom_tiplab(col = color.vector,
              size = 2)  +
  geom_text2(aes(subset=!isTip, label=node))

pdf("dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_tree_original.pdf",
    height = 40,
    width = 16)
mcrA.tree.object
dev.off()


