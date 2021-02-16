#### code/2019_analysis_assembly/sulfur/dsrA/dsr_tree.R ####
# Benjamin D. Peterson




#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ape)
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)



#### Read in tree ####

dsrA.tree.unrooted <- read.newick(file = "dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/dsrA_phylogeny_masked.tree")
dsrA.tree.unrooted$tip.label[grep("ˆ", dsrA.tree.unrooted$tip.label)] <- gsub("ˆ",
                                                                              "u",
                                                                              dsrA.tree.unrooted$tip.label[grep("ˆ", dsrA.tree.unrooted$tip.label)])


#### Visualize unrooted tree ####

pdf("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/original_dsrA_tree.pdf",
    width = 12,
    height = 120)
ggtree(dsrA.tree.unrooted) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()





#### Root tree ####

dsrA.tree <- root(dsrA.tree.unrooted,
                  node = 911)





#### Visualize rooted tree ####

pdf("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/original_dsrA_rooted_tree.pdf",
    width = 12,
    height = 120)
ggtree(dsrA.tree) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()





#### Get list of dsrA genes ####

dsrA.list <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/dsrA_derep_list.txt")






#### Get list of rdsr genes ####

# First isolate those tip labels
endpoints <- c("Nitrospirae_bacterium_RBG_19FT_COMBO_42_15",
               "RBG_16_scaffold_151951")
mrca.of.interest <- getMRCA(dsrA.tree,
                            endpoints)
rdsr.tree <- tree_subset(dsrA.tree,
                         mrca.of.interest,
                         levels_back = 0)
ggtree(rdsr.tree) +
  geom_tiplab(size = 2)

rdsr.list <- rdsr.tree$tip.label[rdsr.tree$tip.label %in% dsrA.list]
writeLines(rdsr.list,
           "dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/dsrA_rev_list.txt")




#### Find reductive dsrA genes ####

red.dsrA <- dsrA.list[which(!(dsrA.list %in% rdsr.list))]




#### Read out list of reductive dsrA genes ####

writeLines(red.dsrA,
           "dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/dsrA_red_list.txt")













#### Analysis of reductive dsrA ####

rm(list = ls())


#### Read in tree ####

dsrA.tree.unrooted <- read.newick(file = "dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/reductive_dsrA_phylogeny/RAxML_bipartitions.dsrA")
dsrA.tree <- midpoint(dsrA.tree.unrooted)


#### Get list of dsrA genes ####

dsrA.list <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/dsrA_derep_list.txt")
reductive.dsrA <- readLines("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA/dsrA_red_list.txt")


#### Read in metadata ####

dsrA.metadata <- read.table("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/reductive_dsrA_phylogeny/ref_dsrA_metadata.tsv",
                            sep = '\t',
                            stringsAsFactors = FALSE)
names(dsrA.metadata) <- c("accessionID", "organism", "taxID")
dsrA.taxonomy <- read.csv("dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/reductive_dsrA_phylogeny/ref_dsrA_taxInfo.tsv",
                          header = FALSE,
                          stringsAsFactors = FALSE)
names(dsrA.taxonomy) <- c("taxID", "division", "genus", "species")

dsr.tax.metadata <- full_join(dsrA.taxonomy,
                              dsrA.metadata)

dsrA.metadata.vector <- paste(dsr.tax.metadata$division, ": ",
                              dsr.tax.metadata$genus, " ",
                              dsr.tax.metadata$species,
                              " (", dsr.tax.metadata$accessionID, ")",
                              sep = "")
names(dsrA.metadata.vector) <- dsrA.metadata$accessionID


#### Get indices for dsrA groups ####

dsrA.EG.index <- which(dsrA.tree$tip.label %in% dsrA.list)
dsrA.ref.index <- grep("\\.", dsrA.tree$tip.label)


#### Get indices for dsrA sequences for which we have reference labels ####

dsr.ref.index.renaming <- which(dsrA.tree$tip.label %in% names(dsrA.metadata.vector))
dsrA.tree$tip.label[dsr.ref.index.renaming] <- dsrA.metadata.vector[dsrA.tree$tip.label[dsr.ref.index.renaming]]



#### Read in and filter coverage values ####

depth.data <- read.csv("dataEdited/2019_analysis_assembly/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE) %>%
  gather(key = metagenomeID,
         value = coverage,
         -1) %>%
  filter(grepl("KMBP005", metagenomeID)) %>%
  spread(key = metagenomeID,
         value = coverage)
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

MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
renaming.vector <- paste(MG.metadata$siteID,
                         MG.metadata$metagenomeID,
                         sep = "_")
names(renaming.vector) <- MG.metadata$metagenomeID






##### Rename MGs with sample info ####

colnames(depth.data.dsrA) <- renaming.vector[colnames(depth.data.dsrA)]
depth.data.dsrA <- depth.data.dsrA %>%
  select(`2A-N_KMBP005A`, `2A-N_KMBP005G`,
         `2A-A_KMBP005B`, `2A-A_KMBP005H`,
         `3A-O_KMBP005F`, `3A-O_KMBP005L`,
         `3A-N_KMBP005E`, `3A-N_KMBP005K`,
         `3A-F_KMBP005D`, `3A-F_KMBP005J`,
         `LOX8_KMBP005C`, `LOX8_KMBP005I`)


#### Generate tree ####

dsrA.tree.object <- ggtree(dsrA.tree) +
  geom_tiplab(size = 2,
              col = color.vector) +
  geom_treescale(x = 0.2, y = 60) +
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
pdf("results/2019_analysis_assembly/sulfur/dsrA_reductive_tree_abund.pdf",
    height = 15,
    width = 8)
gheatmap(dsrA.tree.object,
         depth.data.dsrA, offset=1, width=0.5, font.size=2, 
         colnames_angle=-90, hjust=0)  +
  scale_fill_viridis_c(option="D", name="continuous\nvalue")
dev.off()



#### Make table with cluster IDs ####
cluster1.tree <- tree_subset(dsrA.tree,
                             node = getMRCA(dsrA.tree, c("Sed992Mega19_000000657867_1", "Sed994Meta19_000000022788_1")),
                             levels_back = 0)
cluster1.seqs <- grep("Sed99",
                      cluster1.tree$tip.label,
                      value = TRUE)
cluster1.df <- data.frame(seqID = cluster1.seqs,
                          clusterID = rep("cluster1", length(cluster1.seqs)))

cluster2.tree <- tree_subset(dsrA.tree,
                             node = getMRCA(dsrA.tree, c("Sed996Meta19_000000051661_2", "Sed991Mega19_000001832978_5")),
                             levels_back = 0)
cluster2.seqs <- grep("Sed99",
                      cluster2.tree$tip.label,
                      value = TRUE)
cluster2.df <- data.frame(seqID = cluster2.seqs,
                          clusterID = rep("cluster2", length(cluster2.seqs)))

cluster3.df <- data.frame(seqID = c("Sed992Mega19_000001926553_6", "Sed992Mega19_000000328278_8"),
                          clusterID = "cluster3")

cluster4.tree <- tree_subset(dsrA.tree,
                             node = getMRCA(dsrA.tree, c("Sed994Meta19_000000032633_1", "Sed994Meta19_000000033050_2")),
                             levels_back = 0)
cluster4.seqs <- grep("Sed99",
                      cluster4.tree$tip.label,
                      value = TRUE)
cluster4.df <- data.frame(seqID = cluster4.seqs,
                          clusterID = rep("cluster4", length(cluster4.seqs)))

cluster.df <- rbind(cluster1.df, cluster2.df, cluster3.df, cluster4.df)



#### Write out cluster dataframe ####
write.csv(cluster.df,
          "dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/reductive_dsrA_phylogeny/dsrA_cluster.csv",
          row.names = FALSE)
