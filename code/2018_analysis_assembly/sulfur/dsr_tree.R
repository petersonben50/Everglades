#### code/2018_analysis_assembly/sulfur/dsr/dsr_tree.R ####
# Benjamin D. Peterson




#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ape)
library(ggtree)
library(phangorn)
library(tidyverse)
library(treeio)



#### Read in tree ####

dsrA.tree.unrooted <- read.newick(file = "dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny_red_vs_reverse/RAxML_bipartitions.dsrA")
dsrA.tree.unrooted$tip.label[grep("ˆ", dsrA.tree.unrooted$tip.label)] <- gsub("ˆ",
                                                                              "u",
                                                                              dsrA.tree.unrooted$tip.label[grep("ˆ", dsrA.tree.unrooted$tip.label)])


#### Visualize unrooted tree ####

pdf("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny_red_vs_reverse/original_dsr_tree.pdf",
    width = 12,
    height = 120)
ggtree(dsrA.tree.unrooted) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()






#### Root tree ####

dsrA.tree <- root(dsrA.tree.unrooted,
                  node = 577)





#### Visualize rooted tree ####

pdf("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny_red_vs_reverse/original_dsr_tree_rooted.pdf",
    width = 12,
    height = 120)
ggtree(dsrA.tree) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()





#### Get list of dsrA genes ####

dsrA.list <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny_red_vs_reverse/dsrA_derep_list.txt")






#### Get list of rdsr genes ####

# First isolate those tip labels
endpoints <- c("Nitrospirae_bacterium_RBG_19FT_COMBO_42_15",
               "RBG_16_scaffold_151951")
mrca <- getMRCA(dsrA.tree,
                endpoints)
rdsr.tree <- tree_subset(dsrA.tree,
                         mrca)
ggtree(rdsr.tree) +
  geom_tiplab(size = 2)

rdsr.list <- rdsr.tree$tip.label[rdsr.tree$tip.label %in% dsrA.list]




#### Find reductive dsrA genes ####

red.dsrA <- dsrA.list[which(!(dsrA.list %in% rdsr.list))]




#### Read out list of reductive dsrA genes ####

writeLines(red.dsrA,
           "dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny_red_vs_reverse/dsrA_red_list.txt")













#### Analysis of reductive dsrA ####

rm(list = ls())


#### Read in tree ####

dsrA.tree.unrooted <- read.newick(file = "dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_reductive_phylogeny/RAxML_bipartitions.dsrA")
dsrA.tree <- midpoint(dsrA.tree.unrooted)


#### Get list of dsrA genes ####

dsrA.list <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny_red_vs_reverse/dsrA_derep_list.txt")
reductive.dsrA <- readLines("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny_red_vs_reverse/dsrA_red_list.txt")


#### Read in metadata ####

dsrA.metadata <- read.table("dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_reductive_phylogeny/ref_dsrA_metadata.tsv",
                            sep = '\t',
                            stringsAsFactors = FALSE)
dsrA.metadata.vector <- dsrA.metadata$V2
names(dsrA.metadata.vector) <- dsrA.metadata$V1


#### Get indices for dsrA groups ####

dsrA.EG.index <- which(dsrA.tree$tip.label %in% dsrA.list)
dsrA.ref.index <- grep("\\|", dsrA.tree$tip.label)


#### Fix up names of references ####

dsrA.tree$tip.label[dsrA.ref.index] <- dsrA.tree$tip.label[dsrA.ref.index] %>%
  strsplit("\\|") %>%
  sapply("[", 2)


#### Get indices for dsrA sequences for which we have reference labels ####

dsr.ref.index.renaming <- which(dsrA.tree$tip.label %in% names(dsrA.metadata.vector))
dsrA.tree$tip.label[dsr.ref.index.renaming] <- dsrA.metadata.vector[dsrA.tree$tip.label[dsr.ref.index.renaming]]



#### Read in and filter coverage values ####

depth.data <- read.csv("dataEdited/2018_analysis_assembly/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE)
dsrA.scaffolds <- paste(dsrA.list %>% strsplit("_") %>% sapply("[", 1),
                        dsrA.list %>% strsplit("_") %>% sapply("[", 2),
                        sep = "_")
names(dsrA.list) <- dsrA.scaffolds
depth.data.dsrA <- depth.data %>%
  filter(scaffoldID %in% dsrA.scaffolds)
rownames(depth.data.dsrA) <- dsrA.list[depth.data.dsrA$scaffoldID]
depth.data.dsrA <- depth.data.dsrA %>% select(-scaffoldID)


#### Generate a color vector ####

color.vector <- rep("black", length(dsrA.tree$tip.label))
color.vector[dsrA.EG.index] <- "red"







##### Read in MG metadata ####

MG.metadata <- read_xlsx("metadata/metagenomes/2018_MGs.xlsx")
renaming.vector <- MG.metadata$siteID
names(renaming.vector) <- MG.metadata$metagenomeID






##### Rename MGs with sample info ####

colnames(depth.data.dsrA) <- renaming.vector[colnames(depth.data.dsrA)]




#### Generate tree ####

dsrA.tree.object <- ggtree(dsrA.tree) +
  geom_tiplab(size = 2,
              col = color.vector) +
  geom_treescale(x = 0.2, y = 60) +
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
png("results/2018_analysis_assembly/metabolicProteins/sulfur/dsrA_reductive_tree_abund.png",
    height = 10,
    width = 8,
    units = "in",
    res = 200)
gheatmap(dsrA.tree.object,
         depth.data.dsrA, offset=1, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0)  +
  scale_fill_viridis_c(option="D", name="continuous\nvalue")
dev.off()

