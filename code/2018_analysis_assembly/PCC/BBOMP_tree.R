#### code/2019_analysis_assembly/PCC/PCC_abundance.R ####
# Benjamin D. Peterson




#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ape)
library(ggtree)
library(patchwork)
library(phangorn)
library(readxl)
library(rhmmer)
library(tidyverse)
library(treeio)






#### Read in color vector ####

CBF.vector <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")





#### Read in reference data ####

ref.key <- read_xlsx("/Users/benjaminpeterson/Box/references/proteins/EET/PCC/pcc_omp_key.xlsx")


rename.vector <- ref.key$organismName
rename.vector[which(!is.na(ref.key$geneName))] <- paste(rename.vector[which(!is.na(ref.key$geneName))],
                                                        " - ",
                                                        ref.key$geneName[which(!is.na(ref.key$geneName))],
                                                        sep = "")
names(rename.vector) <- ref.key$accessionID






##### Read in depth info ####
# 
# depth.data <- read.csv("dataEdited/2018_analysis_assembly/metabolicProteins/PCC/BBOMP_phylogeny/PCC_depth_clean.csv",
#                        stringsAsFactors = FALSE) %>%
#   mutate(total = KMBP004A + KMBP004B + KMBP004C +
#            KMBP004D + KMBP004E + KMBP004F) %>%
#   select(scaffoldID, total)
# 




#### Read in tree file ####


BBOMP.tree.unrooted <- read.newick("dataEdited/2018_analysis_assembly/metabolicProteins/PCC/BBOMP_phylogeny/pcc_omp.tree")
BBOMP.tree <- midpoint(BBOMP.tree.unrooted,
                       node.labels = "support")



#### Rename references ####

ref.index <- which(BBOMP.tree$tip.label %in% names(rename.vector))
BBOMP.tree$tip.label[ref.index] <- rename.vector[BBOMP.tree$tip.label[ref.index]]







#### Color vector ####

tree.color.vector <- rep(CBF.vector["bluishgreen"], length(BBOMP.tree$tip.label))
tree.color.vector[ref.index] <- CBF.vector["black"]







#### Plot tree ####

BBOMP.tree.viz <- ggtree(BBOMP.tree,
       aes(x = 0,
           xend = 4)) +
  geom_tiplab(col = tree.color.vector)


BBOMP.tree.viz
