#### code/2019_readBasedAnalyses/nonpareil_data.R ####
# Benjamin D. Peterson



#### Always start with a clean slate ####
setwd("~/Documents/research/Everglades/")
library(gridExtra)
library(Nonpareil)
library(readxl)
library(tidyverse)



#### Read in metadata for sediment metagenomes ####
metadata.seds.2019 <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx") %>%
  filter(grepl("KMBP005", metagenomeID)) 



MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")



#### Set color palette ####
cb.palette <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
colors.to.use <- cb.palette[1:6]
names(colors.to.use) <- unique(metadata.seds.2019$siteID)


#### Read in sample list  ####

sample_list <- metadata.seds.2019 %>%
  mutate(fileLocation = paste("dataEdited/2019_readBasedAnalysis/nonpareil/",
                      metagenomeID,
                      "_NPoutput.npo",
                      sep = "")) %>%
  mutate(colors_for_plot = colors.to.use[siteID]) %>%
  mutate(locationID = paste(metagenomeID,
                             " (", siteID, ")",
                             sep = ""))
  select(fileLocation, locationID, colors_for_plot)
attach(sample_list)



#### Generate nonpareil object and graph ####
pdf("results/2019_readBasedAnalysis/nonpareil_plots.pdf",
    width = 12,
    height = 6)
nps <- Nonpareil.set(files = fileLocation,
                     col = colors_for_plot,
                     labels = locationID)
dev.off()
detach(sample_list)



#### Plot diversity ####
np.diversity <- summary(nps)[,"diversity"]
np.diversity.df <- data.frame(locationID = names(np.diversity),
                              diversity = np.diversity) %>%
  mutate(siteID = locationID %>%
           strsplit("\\(") %>% sapply("[", 2) %>%
           strsplit("\\)") %>% sapply("[", 1))

pdf("results/2019_readBasedAnalysis/nonpareil_diversity.pdf",
    width = 6,
    height = 4)
np.diversity.df %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  ggplot(aes(x = siteID,
             y = diversity)) +
  geom_point() +
  theme_classic() +
  ylim(c(20, 25)) +
  ylab("Nonpareil diversity index")
dev.off()

# Run a two-way ANOVA on diversity
two.way.anova <- aov(diversity ~ siteID,
                     data = np.diversity.df)
sink("results/2019_readBasedAnalysis/nonpareil_diversity_anova.txt")
summary(two.way.anova)
sink()
# Not significantly different.



#### Make chart of coverage ####
MGcoverage.vector <- summary(nps)[,"C"]*100
MGcoverage.df <- data.frame(siteID = names(MGcoverage.vector) %>%
                              strsplit("\\(") %>% sapply("[", 2) %>%
                              strsplit("\\)") %>% sapply("[", 1),
                            metagenomeID = names(MGcoverage.vector) %>%
                              strsplit(" \\(") %>% sapply("[", 1),
                            percent_coverage = MGcoverage.vector) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  arrange(siteID)
pdf("results/2019_readBasedAnalysis/nonpareil_coverage.pdf",
    width = 4,
    height = 4)
grid.table(MGcoverage.df,
           rows = NULL,
           cols = colnames(MGcoverage.df))
dev.off()


#### Also save to csv ####
write.csv(MGcoverage.df,
          "results/2019_readBasedAnalysis/nonpareil_coverage.csv",
          row.names = FALSE)
