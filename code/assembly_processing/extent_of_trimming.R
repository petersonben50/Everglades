#### code/assembly_processing/extent_of_trimming.R ####
# Benjamin D. Peterson

# This script will read in the coverage of the
# metagenomes from the Everglades. It will read
# out a vector for normalizing the coverages of
# sequences in those metagenomes.



#### Get your junk out of here ####

rm(list = ls())
setwd("Documents/research/Everglades/")
library(tidyverse)



#### Read in total coverage data pre-trimming ####


coverage.data.pre.trim <- read.table("dataEdited/metagenomes/reports/metagenome_coverage_pre_trimming.tsv",
                                     sep = '\t',
                                     header = TRUE,
                                     stringsAsFactors = FALSE) %>%
  mutate(total_coverage_pre = R1 + R2) %>%
  select(metagenomeID, total_coverage_pre)





#### Read in total coverage data post-trimming ####

coverage.data.post.trim <- read.table("dataEdited/metagenomes/reports/metagenome_coverage.tsv",
                                      sep = '\t',
                                      header = TRUE,
                                      stringsAsFactors = FALSE) %>%
  mutate(total_coverage_post = R1 + R2 + single + merged) %>%
  select(metagenomeID, total_coverage_post)




#### Check out difference in coverage ####

difference.in.coverage <- full_join(coverage.data.pre.trim,
                                    coverage.data.post.trim) %>%
  mutate(per_change_coverage = (total_coverage_post - total_coverage_pre) / total_coverage_pre * 100)


rm(list = ls(pattern = "coverage"))




#### Read in read counts, pre-trimming ####

read.counts.pre <- read.table("dataEdited/metagenomes/reports/metagenome_read_count_pre_trimming.tsv",
                              sep = '\t',
                              header = TRUE,
                              stringsAsFactors = FALSE)







#### Read in read counts, post-trimming ####

read.counts <- read.table("dataEdited/metagenomes/reports/metagenome_read_count.tsv",
                          sep = '\t',
                          header = TRUE,
                          stringsAsFactors = FALSE)
