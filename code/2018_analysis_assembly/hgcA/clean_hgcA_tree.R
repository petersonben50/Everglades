#### code/2018_analysis_assembly/hgcA/clean_hgcA_tree.R ####
# Benjamin D. Peterson

# This script will look at our phylogenetic tree 
# iterations and ultimately clean up the final one.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ggtree)
library(phangorn)
library(tidyverse)
library(treeio)


#### List of hgcA sequences ####

hgcA.list <- readLines("dataEdited/2018_analysis_assembly/hgcA/hgcA.txt")



#### CB-friendly color vector ####

colorblind.color.vector <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Vectors for confirmed methylators ####

hgcA.references.confirmed.vector <- readRDS("/Users/benjaminpeterson/Box/references/hgcA/hgcA_confirmed_ref_renaming_vector.rds")
ref.confirmed.colors <- readRDS("/Users/benjaminpeterson/Box/references/hgcA/hgcA_confirmed_ref_color_vector.rds")


#### Check out FastTree with just McDaniel 2020 references ####

# Read in tree
tree.name <- "dataEdited/2018_analysis_assembly/hgcA/phylogeny/rough_hgcA.tree"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


# Root tree
hgcA.tree <- midpoint(hgcA.tree.unrooted,
                      node.labels = "support")


# Read in reference data
hgcA.references.McD <- read.csv("/Users/benjaminpeterson/Box/references/hgcA/McDaniel_2020/McDaniel_2020_hgcA_ref_metadata.csv",
                            stringsAsFactors = FALSE) %>%
  select(genome_name, name) %>%
  mutate(tree_name = paste(strsplit(name, " ") %>%
                             sapply("[", 1),
                           " ",
                           strsplit(name, " ") %>%
                             sapply("[", 2),
                           " (",
                           genome_name,
                           ")",
                           sep = ""))
hgcA.references.McD.vector <- hgcA.references.McD$tree_name
names(hgcA.references.McD.vector) <- hgcA.references.McD$genome_name


# Get indices
reference.indices <- which(hgcA.tree$tip.label %in% names(hgcA.references.McD.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)


# Change names:
hgcA.tree$tip.label[reference.indices] <- hgcA.references.McD.vector[hgcA.tree$tip.label[reference.indices]]


# Set color vector
color.vector <- rep("grey", length(hgcA.tree$tip.label))
color.vector[this.study.indices] <- "red"
color.vector[reference.indices] <- "black"

# View tree
ggtree(hgcA.tree, aes(x = 0,
                      xend = 2)) + 
  geom_tiplab(size=2.5, align = TRUE, col = color.vector) + 
  geom_nodelab(aes(x = branch), vjust = -.3)






#### Check out FastTree with McDaniel 2020 and NCBI non-redundant references ####

# Read in tree
tree.name <- "dataEdited/2018_analysis_assembly/hgcA/phylogeny/rough_hgcA_2.tree"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


# Root tree
hgcA.tree <- midpoint(hgcA.tree.unrooted,
                      node.labels = "support")


# Read in reference data for nr database
hgcA.references.nr <- read.table("references/hgcA/hgcA_ref_nr_taxonomy.tsv",
                                  stringsAsFactors = FALSE,
                                  sep = '\t') %>%
  rename(accessionID = V1,
         name = V2) %>%
  mutate(tree_name = paste(strsplit(name, " ") %>%
                             sapply("[", 1),
                           " ",
                           strsplit(name, " ") %>%
                             sapply("[", 2),
                           " (",
                           accessionID,
                           ")",
                           sep = "")) %>% 
  select(accessionID, tree_name)

hgcA.references.nr.vector <- hgcA.references.nr$tree_name
names(hgcA.references.nr.vector) <- hgcA.references.nr$accessionID



# Get indices
reference.indices.McD <- which(hgcA.tree$tip.label %in% names(hgcA.references.McD.vector))
reference.indices.nr <- which(hgcA.tree$tip.label %in% names(hgcA.references.nr.vector))
reference.indices.confirmed <- which(hgcA.tree$tip.label %in% names(hgcA.references.confirmed.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)


# Change names:
hgcA.tree$tip.label[reference.indices.McD] <- hgcA.references.McD.vector[hgcA.tree$tip.label[reference.indices.McD]]
hgcA.tree$tip.label[reference.indices.nr] <- hgcA.references.nr.vector[hgcA.tree$tip.label[reference.indices.nr]]
hgcA.tree$tip.label[reference.indices.confirmed] <- hgcA.references.confirmed.vector[hgcA.tree$tip.label[reference.indices.confirmed]]


# Set color vector
color.vector <- rep("grey", length(hgcA.tree$tip.label))
color.vector[c(reference.indices.McD, reference.indices.nr)] <- colorblind.color.vector["black"]
color.vector[reference.indices.confirmed] <- colorblind.color.vector[ref.confirmed.colors[hgcA.tree$tip.label[reference.indices.confirmed]]]
color.vector[this.study.indices] <- "red"

# View tree
pdf("results/2018_assembly_analysis/hgcA/phylogeny/rough_tree_2.pdf",
    width = 5,
    height = 15)
ggtree(hgcA.tree, aes(x = 0,
                      xend = 4)) + 
  geom_tiplab(size=2.5, align = TRUE, col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 1)
dev.off()










#### Check out FastTree of final alignment ####

# Read in tree
tree.name <- "dataEdited/2018_analysis_assembly/hgcA/phylogeny/rough_hgcA_final.tree"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


# Root tree
hgcA.tree <- midpoint(hgcA.tree.unrooted,
                      node.labels = "support")

# Get indices
reference.indices.McD <- which(hgcA.tree$tip.label %in% names(hgcA.references.McD.vector))
reference.indices.nr <- which(hgcA.tree$tip.label %in% names(hgcA.references.nr.vector))
reference.indices.confirmed <- which(hgcA.tree$tip.label %in% names(hgcA.references.confirmed.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)

# Change names:
hgcA.tree$tip.label[reference.indices.McD] <- hgcA.references.McD.vector[hgcA.tree$tip.label[reference.indices.McD]]
hgcA.tree$tip.label[reference.indices.nr] <- hgcA.references.nr.vector[hgcA.tree$tip.label[reference.indices.nr]]
hgcA.tree$tip.label[reference.indices.confirmed] <- hgcA.references.confirmed.vector[hgcA.tree$tip.label[reference.indices.confirmed]]

# Set color vector
color.vector <- rep("grey", length(hgcA.tree$tip.label))
color.vector[c(reference.indices.McD, reference.indices.nr)] <- colorblind.color.vector["black"]
color.vector[reference.indices.confirmed] <- colorblind.color.vector[ref.confirmed.colors[hgcA.tree$tip.label[reference.indices.confirmed]]]
color.vector[this.study.indices] <- "red"

# Visualize tree
ggtree(hgcA.tree, aes(x = 0,
                      xend = 4)) + 
  geom_tiplab(size=2.5, align = TRUE, col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 1)




#### Clean up RAxML tree from final alignment ####

# Read in tree
tree.name <- "dataEdited/2018_analysis_assembly/hgcA/phylogeny/RAxML_bipartitions.hgcA"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


# Root tree
hgcA.tree <- midpoint(hgcA.tree.unrooted,
                      node.labels = "support")

# Get indices
reference.indices.McD <- which(hgcA.tree$tip.label %in% names(hgcA.references.McD.vector))
reference.indices.nr <- which(hgcA.tree$tip.label %in% names(hgcA.references.nr.vector))
reference.indices.confirmed <- which(hgcA.tree$tip.label %in% names(hgcA.references.confirmed.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)

# Change names:
hgcA.tree$tip.label[reference.indices.McD] <- hgcA.references.McD.vector[hgcA.tree$tip.label[reference.indices.McD]]
hgcA.tree$tip.label[reference.indices.nr] <- hgcA.references.nr.vector[hgcA.tree$tip.label[reference.indices.nr]]
hgcA.tree$tip.label[reference.indices.confirmed] <- hgcA.references.confirmed.vector[hgcA.tree$tip.label[reference.indices.confirmed]]

# Set color vector
color.vector <- rep("grey", length(hgcA.tree$tip.label))
color.vector[c(reference.indices.McD, reference.indices.nr)] <- colorblind.color.vector["black"]
color.vector[reference.indices.confirmed] <- colorblind.color.vector[ref.confirmed.colors[hgcA.tree$tip.label[reference.indices.confirmed]]]
color.vector[this.study.indices] <- "red"

# Visualize tree
ggtree(hgcA.tree, aes(x = 0,
                      xend = 3)) + 
  geom_tiplab(size=2.5, align = TRUE, col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)



#### Save out tree ####
saveRDS(hgcA.tree,
        "dataEdited/2018_analysis_assembly/hgcA/phylogeny/hgcA_clean_tree.rds")


#### Save out color vector ####
saveRDS(color.vector,
        "dataEdited/2018_analysis_assembly/hgcA/phylogeny/hgcA_clean_tree_color_vector.rds")
