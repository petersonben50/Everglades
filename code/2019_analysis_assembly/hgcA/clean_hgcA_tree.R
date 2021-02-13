#### code/2019_analysis_assembly/hgcA/clean_hgcA_tree.R ####
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


#### List of hgcA sequences ####

hgcA.list <- readLines("dataEdited/2019_analysis_assembly/hgcA/hgcA.txt")



#### CB-friendly color vector ####

colorblind.color.vector <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Vectors for confirmed methylators ####

hgcA.references.confirmed.vector <- readRDS("/Users/benjaminpeterson/Box/references/hgcA/hgcA_confirmed_ref_renaming_vector.rds")
ref.confirmed.colors <- readRDS("/Users/benjaminpeterson/Box/references/hgcA/hgcA_confirmed_ref_color_vector.rds")


#### Check out FastTree of final alignment ####

# Read in tree
tree.name <- "dataEdited/2019_analysis_assembly/hgcA/phylogeny/rough_hgcA.tree"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


# Root tree
hgcA.tree <- midpoint(hgcA.tree.unrooted,
                      node.labels = "support")

# Get indices
reference.indices.confirmed <- which(hgcA.tree$tip.label %in% names(hgcA.references.confirmed.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)

# Change names:
hgcA.tree$tip.label[reference.indices.confirmed] <- hgcA.references.confirmed.vector[hgcA.tree$tip.label[reference.indices.confirmed]]

# Set color vector
color.vector <- rep("grey", length(hgcA.tree$tip.label))
color.vector[reference.indices.confirmed] <- colorblind.color.vector[ref.confirmed.colors[hgcA.tree$tip.label[reference.indices.confirmed]]]
color.vector[this.study.indices] <- "red"

# Visualize tree
pdf("results/2019_analysis_assembly/hgcA_tree_FastTree.pdf",
    height = 120,
    width = 8)
ggtree(hgcA.tree) + 
  geom_tiplab(size=2.5, align = TRUE, col = color.vector) +
  geom_text2(aes(subset=!isTip, label=node),
             size = 1.5)
dev.off()



#### Pull out sequence names to remove ####

nodes.to.remove <- c(1671, 1766, 1599, 1619, 1626,
                     1628, 1472, 1415, 1388, 1376,
                     1537, 1372, 1545, 1157, 1267,
                     1286, 1137, 1295, 1845, 1829,
                     1953, 1984, 1999, 1941, 1929,
                     1919, 1906, 1897, 2052, 2045,
                     2039, 2030, 2083, 2087, 2097,
                     2090, 2102, 2112, 2124, 2126,
                     2128)
seqs.to.remove <- vector()
for (node.of.interest in nodes.to.remove) {
  tree.subset.to.remove <- tree_subset(hgcA.tree,
                                       node = node.of.interest,
                                       levels_back = 0)
  seqs.to.remove <- c(seqs.to.remove, tree.subset.to.remove$tip.label)
}
seqs.to.remove <- grep("_ISO_",
                       seqs.to.remove,
                       invert = TRUE,
                       value = TRUE)

fileConn <- file("dataEdited/2019_analysis_assembly/hgcA/phylogeny/seqs_to_remove.txt")
writeLines(seqs.to.remove, fileConn)
close(fileConn)

# Clean up
rm(list = ls())




#### Clean up RAxML tree from final alignment ####


#### Read in needed files ####

hgcA.list <- readLines("dataEdited/2019_analysis_assembly/hgcA/hgcA.txt")
colorblind.color.vector <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
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
tree.name <- "dataEdited/2019_analysis_assembly/hgcA/phylogeny/RAxML_bipartitions.hgcA"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)
#### Fix tip labels ####
hgcA.tree.unrooted$tip.label <- paste(hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 1),
                                      hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 2),
                                      hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 3),
                                      sep = "_")

# Root tree
hgcA.tree <- root(phy = hgcA.tree.unrooted,
                  outgroup = c("paralog_Thermosulfurimonas_dismutans",
                               "paralog_Candidatus_Omnitrophica",
                               "Sed993Meta19_000000003033_5",
                               "Sed993Meta19_000000001064_12",
                               "Sed993Meta19_000000027668_3",
                               "Sed993Mega19_000001097352_5",
                               "Sed993Meta19_000000035661_4"),
                  edgelabel = TRUE)


# Get indices
reference.indices.mate <- which(hgcA.tree$tip.label %in% names(MATE.renaming.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)

# Change names:
hgcA.tree$tip.label[reference.indices.mate] <- MATE.renaming.vector[hgcA.tree$tip.label[reference.indices.mate]]


# Set color vector
color.vector <- rep("grey", length(hgcA.tree$tip.label))
color.vector[reference.indices.mate] <- colorblind.color.vector["black"]
color.vector[this.study.indices] <- colorblind.color.vector["vermillion"]

# Visualize tree
pdf("results/2019_analysis_assembly/hgcA_tree_RAxML_rooted.pdf",
    height = 60,
    width = 10)
ggtree(hgcA.tree, aes(x = 0,
                      xend = 4.5)) + 
  geom_tiplab(size=2.5, align = TRUE, col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
dev.off()




#### Save out tree ####
saveRDS(hgcA.tree,
        "dataEdited/2019_analysis_assembly/hgcA/phylogeny/hgcA_clean_tree.rds")


#### Save out color vector ####
saveRDS(color.vector,
        "dataEdited/2019_analysis_assembly/hgcA/phylogeny/hgcA_clean_tree_color_vector.rds")
