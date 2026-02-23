# Beta Diveristy Analysis 
library(phyloseq)
library(tidyverse)
library(vegan)

# Loading Object
ps = readRDS("Datasets/phyloseq_taxonomy.rds")

gastric_metadata <- read.delim("Datasets/gastric_cancer_metadata.tsv", row.names = 1)


# Rarefy the Data

psrare = ps %>% rarefy_even_depth(sample.size = 8000, rngseed = 421)

# Distance between each sample pair calculation
ps_bray = phyloseq::distance(psrare, method = "bray")
as.matrix(ps_bray)

# MDS scaling 
set.seed(421)

mds = metaMDS(ps_bray)

# Extracting Data

mds_data = mds$points %>% as.data.frame %>%
  merge(sample_data(psrare), by='row.names', sort=F)
head(mds_data)


# Plotting ~ Separating Each Group and Comparing Gastric Microbiome Based on 
# Biological Sex

beta_diversity_plot = mds_data %>%
  ggplot(aes(MDS1,MDS2,color = Gender)) +
  geom_point() +
  stat_ellipse() +
  facet_wrap(~Group) +
  theme_classic(base_size=18) +
  labs(
    title = "Bray–Curtis Beta Diversity by Sex Within Each Gastric Disease Stage",
    x = "MDS1",
    y = "MDS2",
    color = "Biological Sex") +
  theme_bw()

beta_diversity_plot

# Note: Need to check the stratification of the overall data to determine which
# variable is pulling a lot weight (i.e., driving overall microbiome variation)

# Saving the Plot!
ggsave("Results/Plots/beta_diversity_bray-curtis.png",
       beta_diversity_plot,
       width= 10,
       height = 10)



#### PERMANOVA Analysis ####

# Select GC gastric disease stage to analyze
disease_stage <- c("Gastric cancer (GC)") # Repeat this with the four other stages!

# Kepp ONLY samples from that stage
samples_to_keep = gastric_metadata %>% filter(Group %in% disease_stage) %>% 
  rownames()

# Identify samples within the selected gastric disease stage that are present in 
# the Bray-Curtis distance matrix

dm_ids <- rownames(as.matrix(ps_bray))

samples_to_keep <- gastric_metadata %>%
  filter(Group == disease_stage) %>%
  rownames() %>%
  intersect(dm_ids)

ps_bray_sub <- as.matrix(ps_bray)[samples_to_keep, samples_to_keep] %>% 
  as.dist()

# Obtaining sample IDs that are present in the Bray-Curtis distance matrix
keep_ids <-  rownames(as.matrix(ps_bray_sub))

# Subsetting the metadata so that only those samples are present for analysis in 
# the same order
metadata_sub <- gastric_metadata [keep_ids, , drop = FALSE]


# Run PERMANOVA Comparing Biological Sex within GC Stage with the subsetted metadata

stats_stage <- adonis2(ps_bray_sub ~ Gender, data = metadata_sub)

stats_stage

# To Do: Continue on the PERMANOVA with other stages 





