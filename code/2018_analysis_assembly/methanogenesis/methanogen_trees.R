#### code/2018_analysis_assembly/methanogenesis/methanogen_trees.R ####
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

mcrA.tree.unrooted <- read.newick("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_for_phylogeny.tree")
mcrA.tree <- midpoint(mcrA.tree.unrooted)


#### Read in mcrA list from this study ####

mcrA.list <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA/mcrA_derep_list.txt")


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

pdf("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_tree_original.pdf",
    height = 40,
    width = 16)
mcrA.tree.object
dev.off()



#### Nodes to remove ####

nodes.to.remove <- c(263, 260, 250, 248, 269, 242,
                     276, 279, 336, 358, 402, 401,
                     382, 394, 403, 432)
tips.to.remove <- NULL
for (node.of.interest in nodes.to.remove) {
  small.tree <- tree_subset(mcrA.tree,
                            node = node.of.interest,
                            levels_back = 0)
  tips.to.remove <<- c(small.tree$tip.label, tips.to.remove)
}
additional.seqs.to.remove <- c("SRR1367222_1", "SRR634686_1", "SRR1010317_1", "mgm4537093_2",
                               "ClassI−15_Methanothermus_fervidus_WP_013413861.1")
tips.to.remove <- c(tips.to.remove, additional.seqs.to.remove)


fileConn<-file("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/seqs_to_remove.txt")
writeLines(tips.to.remove, fileConn)
close(fileConn)

rm(list = ls())



















#### Read in mcrA tree ####

mcrA.tree.unrooted <- read.newick("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_for_phylogeny_2.tree")
mcrA.tree <- midpoint(mcrA.tree.unrooted)


#### Read in mcrA list from this study ####

mcrA.list <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA/mcrA_derep_list.txt")


#### Make indices ####

my.seq.index <- which(mcrA.tree$tip.label %in% mcrA.list)


#### Make color vector ####

color.vector <- rep("black", length(mcrA.tree$tip.label))
color.vector[my.seq.index] <- "red"

#### Read in and filter coverage values ####

depth.data <- read.csv("dataEdited/2018_analysis_assembly/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE)
mcrA.scaffolds <- paste(mcrA.list %>% strsplit("_") %>% sapply("[", 1),
                        mcrA.list %>% strsplit("_") %>% sapply("[", 2),
                        sep = "_")
names(mcrA.list) <- mcrA.scaffolds
depth.data.mcrA <- depth.data %>%
  filter(scaffoldID %in% mcrA.scaffolds)
rownames(depth.data.mcrA) <- mcrA.list[depth.data.mcrA$scaffoldID]
depth.data.mcrA <- depth.data.mcrA %>% select(-scaffoldID)


#### Make mcrA tree ####


mcrA.tree.object <- ggtree(mcrA.tree) +
  geom_tiplab(col = color.vector,
              size = 2)  +
  geom_nodelab()

mcrA.tree.object


rm(list = ls(pattern = "mcrA.tree"),
   color.vector, my.seq.index)




#### Read in RAxML tree ####

mcrA.tree.unrooted <- read.newick("dataEdited/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/take_2/RAxML_bipartitions.mcrA")
outgroup.ID <- "OFV68676_1_S_caldarius_mcrA1"
mcrA.tree <- root(mcrA.tree.unrooted,
                  outgroup = outgroup.ID)


#### Subset tree to remove root ####

mcrA.tree <- drop.tip(mcrA.tree, outgroup.ID)



#### Make indices ####

my.seq.index <- which(mcrA.tree$tip.label %in% mcrA.list)


#### Make color vector ####

color.vector <- rep("black", length(mcrA.tree$tip.label))
color.vector[my.seq.index] <- "red"



#### Make final tree object ####
mcrA.tree.object <- ggtree(mcrA.tree) +
  geom_tiplab(col = color.vector,
              size = 2)

png("results/2018_analysis_assembly/metabolicProteins/methanogenesis/mcrA_tree_abund.png",
    res = 200,
    units = "in",
    height = 14,
    width = 8)
gheatmap(mcrA.tree.object,
         depth.data.mcrA, offset=1, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0)  +
  scale_fill_viridis_c(option="D", name="continuous\nvalue")
dev.off()
