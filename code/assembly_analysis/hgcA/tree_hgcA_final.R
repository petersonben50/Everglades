#### code/assembly_analysis/hgcA/tree_hgcA_final.R ####
# Benjamin D. Peterson

# This script will look at our phylogenetic tree 
# iterations and ultimately clean up the final one.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)


#### Clean up RAxML tree from final alignment ####


#### Read in needed files ####

hgcA.list <- readLines("dataEdited/assembly_analysis/hgcA/hgcA.txt")
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
hgcA.references.confirmed.vector <- readRDS("/Users/benjaminpeterson/Box/references/hgcA/hgcA_confirmed_ref_renaming_vector.rds")
ref.confirmed.colors <- readRDS("/Users/benjaminpeterson/Box/references/hgcA/hgcA_confirmed_ref_color_vector.rds")

#### Read in Hg-MATE seq info ####
hg.mate.metadata <- read_xlsx("/Users/benjaminpeterson/Documents/research/Hg_MATE/versions/v1.01142021/Hg-MATE-Db.v1.01142021_catalogue.xlsx") %>%
  rename(group = `microbial group`,
         name = `[Organism Name]_Phylum-Class`) %>%
  mutate(name = gsub("sp._", "sp.", name)) %>%
  mutate(name = paste(name %>% strsplit("_") %>% sapply("[", 1),
                      name %>% strsplit("_") %>% sapply("[", 2),
                      sep = "_"),
         treeID = paste(group,
                        " (", name, "-", MATE_id, ")",
                        sep = ""))
MATE.renaming.vector <- hg.mate.metadata$treeID
names(MATE.renaming.vector) <- hg.mate.metadata$MATE_id


# Read in tree
tree.name <- "dataEdited/assembly_analysis/hgcA/phylogeny/final/RAxML_bipartitions.hgcA"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)

#### Fix tip labels ####
hgcA.tree.unrooted$tip.label <- paste(hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 1),
                                      hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 2),
                                      hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 3),
                                      sep = "_")

#### Look for branch to paralogs ####
pdf("dataEdited/assembly_analysis/hgcA/phylogeny/final/initial_tree.pdf",
    height = 40,
    width = 12)
ggtree(hgcA.tree.unrooted) + 
  geom_tiplab(size=2.5, align = TRUE) +
  geom_text2(aes(subset=!isTip, label=node),
             size = 2.5)
dev.off()
# Node leading to paralogs appears to be 306

# Root tree
hgcA.tree.rooted <- root(phy = hgcA.tree.unrooted,
                         node = 306,
                         edgelabel = TRUE)
ggtree(hgcA.tree.rooted) + 
  geom_tiplab(size=2.5, align = TRUE) +
  geom_text2(aes(subset=!isTip, label=node),
             size = 1.5)
# This looks decent. Those fused proteins are pretty divergent though.

#### Remove root ####
hgcA.tree <- drop.tip(hgcA.tree.rooted,
                      c("Sed993Meta19_000000035661_4",
                        "Sed993Mega19_000001097352_5",
                        "Sed993Meta19_000000027668_3",
                        "Sed993Meta19_000000001064_12",
                        "Sed993Meta19_000000003033_5"))
ggtree(hgcA.tree) + 
  geom_tiplab(size=2.5, align = TRUE) +
  geom_text2(aes(subset=!isTip, label=node),
             size = 1.5)


# Get indices
reference.indices.mate <- which(hgcA.tree$tip.label %in% names(MATE.renaming.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)

# Change names:
hgcA.tree$tip.label[reference.indices.mate] <- MATE.renaming.vector[hgcA.tree$tip.label[reference.indices.mate]]


# Set color vector
color.vector <- rep("grey", length(hgcA.tree$tip.label))
color.vector[reference.indices.mate] <- cb.translator["black"]
color.vector[this.study.indices] <- cb.translator["vermillion"]

# Visualize tree
pdf("results/metagenomes/assembly/hgcA/hgcA_tree_RAxML_final.pdf",
    height = 20,
    width = 10)
ggtree(hgcA.tree,
       aes(x = 0, xend = 4.5)) + 
  geom_tiplab(size=2.5, align = TRUE, colour = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
dev.off()




#### Save out tree ####
saveRDS(hgcA.tree,
        "dataEdited/assembly_analysis/hgcA/phylogeny/final/hgcA_clean_tree.rds")


#### Save out color vector ####
saveRDS(color.vector,
        "dataEdited/assembly_analysis/hgcA/phylogeny/final/hgcA_clean_tree_color_vector.rds")
