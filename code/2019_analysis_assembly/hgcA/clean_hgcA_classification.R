#### code/2019_analysis_assembly/hgcA/clean_hgcA_classification.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(BoSSA)
library(phyloseq)
library(tidyverse)



#### Read in data ####
refpkg_path <- "references/hgcA/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg"
taxonomy <- refpkg(refpkg_path,
                   type = "taxonomy")
sqlite_file <- "dataEdited/2019_analysis_assembly/hgcA/classification/output/Hg_MATE_classify"
jplace_file <- "dataEdited/2019_analysis_assembly/hgcA/classification/output/hgcA_for_classification.jplace"
pplace_object <- read_sqlite(sqlite_file,
                             jplace_file)



#### Generate taxonomy table ####
taxonomy.table <- pplace_to_taxonomy(pplace_object,
                                     taxonomy,
                                     tax_name = TRUE,
                                     rank = c("phylum", "class", "order",
                                              "family", "genus", "species"))
taxonomy.table.df <- as.data.frame(taxonomy.table) %>%
  mutate(seqID = row.names(taxonomy.table)) %>%
  select(seqID, phylum, class, order,
         family, genus, species) %>%
  mutate(phylum = phylum %>%
           strsplit(" <") %>% sapply("[", 1))



#### Save it out ####
write.csv(taxonomy.table.df,
          "dataEdited/2019_analysis_assembly/hgcA/classification/hgcA_taxonomy_table.csv",
          row.names = FALSE)
