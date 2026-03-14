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
  geom_point(size = 3) +
  stat_ellipse(linewidth = 1.1) +
  facet_wrap(~Group, nrow = 2) +
  theme_classic(base_size=18) +
  labs(
    x = "MDS1",
    y = "MDS2",
    color = "Biological Sex") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16.5, face = "bold"),
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    panel.spacing = unit(2.5, "lines")
  ) 

beta_diversity_plot


# Saving the Plot!
ggsave("Results/Plots/beta_diversity_bray-curtis.png",
       beta_diversity_plot,
       width= 10,
       height = 10)

# Note: Need to check the stratification of the overall data to determine which
# variable is pulling a lot weight (i.e., driving overall microbiome variation)

#### Stratification/Overview ####

# Sample IDs in Bray-Curtis distance matrix
bray_ids <- rownames(as.matrix(ps_bray))

# Filtering metadata so it matches the samples used in the distance matrix
gastric_data_bray <- gastric_metadata[bray_ids, ,drop = FALSE]

# Stratification 
stratification <- adonis2(ps_bray ~ Group * Gender, data = gastric_data_bray)
stratification

# To assess whether gastric cancer stage, biological sex, or the interaction
# between these two variables overall variation in the gastric microbiome. 

# The overall PERMANOVA is significant with p = 0.001 and explains ~26.6% variation
# in the microbiome, which indicates that gastric disease stage, biological sex, 
# and/or their interactions have a role in differences observed in gastric
# microbiome beta diversity throughout the dataset. This suggests that these variables
# may play a role in shaping microbial community composition. This overall analysis 
# provides context for the subsequent stage-specific PERMANOVA analyses. 

#### PERMANOVA Analysis ####

# Select GC gastric disease stage to analyze
gc_stage <- c("Gastric cancer (GC)") # Repeat this with the four other stages!

# Keep ONLY samples from that stage
gc_samples_to_keep = gastric_metadata %>% filter(Group %in% gc_stage) %>% 
  rownames()

# Identify samples within the selected gastric disease stage that are present in 
# the Bray-Curtis distance matrix

gc_dm_ids <- rownames(as.matrix(ps_bray))

gc_samples_to_keep <- gastric_metadata %>%
  filter(Group == gc_stage) %>%
  rownames() %>%
  intersect(gc_dm_ids)

gc_bray_sub <- as.matrix(ps_bray)[gc_samples_to_keep, gc_samples_to_keep] %>% 
  as.dist()

# Obtaining sample IDs that are present in the Bray-Curtis distance matrix
gc_keep_ids <-  rownames(as.matrix(gc_bray_sub))

# Subsetting the metadata so that only those samples are present for analysis in 
# the same order
gc_metadata_sub <- gastric_metadata [gc_keep_ids, , drop = FALSE]


# Run PERMANOVA Comparing Biological Sex within GC Stage with the subsetted metadata

stats_gc <- adonis2(gc_bray_sub ~ Gender, data = gc_metadata_sub)

stats_gc

# To Do: Continue on the PERMANOVA with other stages 

## Healthy Control (HC) ##

# Select HC gastric disease stage to analyze
hc_stage <- c("Healthy control (HC)") 

# Keep ONLY samples from that stage
hc_samples_to_keep = gastric_metadata %>% filter(Group %in% hc_stage) %>% 
  rownames()

# Identify samples within the selected gastric disease stage that are present in 
# the Bray-Curtis distance matrix

hc_dm_ids <- rownames(as.matrix(ps_bray))

hc_samples_to_keep <- gastric_metadata %>%
  filter(Group == hc_stage) %>%
  rownames() %>%
  intersect(hc_dm_ids)

hc_bray_sub <- as.matrix(ps_bray)[hc_samples_to_keep, hc_samples_to_keep] %>% 
  as.dist()

# Obtaining sample IDs that are present in the Bray-Curtis distance matrix
hc_keep_ids <-  rownames(as.matrix(hc_bray_sub))

# Subsetting the metadata so that only those samples are present for analysis in 
# the same order
hc_metadata_sub <- gastric_metadata [hc_keep_ids, , drop = FALSE]


# Run PERMANOVA Comparing Biological Sex within HC Stage with the subsetted metadata

stats_hc <- adonis2(hc_bray_sub ~ Gender, data = hc_metadata_sub)

stats_hc

## Intestinal metaplasia (IM) ##

# Select IM gastric disease stage to analyze
im_stage <- c("Intestinal metaplasia (IM）") 

# Keep ONLY samples from that stage
im_samples_to_keep = gastric_metadata %>% filter(Group %in% im_stage) %>% 
  rownames()

# Identify samples within the selected gastric disease stage that are present in 
# the Bray-Curtis distance matrix

im_dm_ids <- rownames(as.matrix(ps_bray))

im_samples_to_keep <- gastric_metadata %>%
  filter(Group == im_stage) %>%
  rownames() %>%
  intersect(im_dm_ids)

im_bray_sub <- as.matrix(ps_bray)[im_samples_to_keep, im_samples_to_keep] %>% 
  as.dist()

# Obtaining sample IDs that are present in the Bray-Curtis distance matrix
im_keep_ids <-  rownames(as.matrix(im_bray_sub))

# Subsetting the metadata so that only those samples are present for analysis in 
# the same order
im_metadata_sub <- gastric_metadata [im_keep_ids, , drop = FALSE]


# Run PERMANOVA Comparing Biological Sex within IM Stage with the subsetted metadata

stats_im <- adonis2(im_bray_sub ~ Gender, data = im_metadata_sub)

stats_im

## Intraepithelial neoplasia (IN) ##

# Select IN gastric disease stage to analyze
in_stage <- c("Intraepithelial neoplasia (IN)") 

# Keep ONLY samples from that stage
in_samples_to_keep = gastric_metadata %>% filter(Group %in% in_stage) %>% 
  rownames()

# Identify samples within the selected gastric disease stage that are present in 
# the Bray-Curtis distance matrix

in_dm_ids <- rownames(as.matrix(ps_bray))

in_samples_to_keep <- gastric_metadata %>%
  filter(Group == in_stage) %>%
  rownames() %>%
  intersect(in_dm_ids)

in_bray_sub <- as.matrix(ps_bray)[in_samples_to_keep, in_samples_to_keep] %>% 
  as.dist()

# Obtaining sample IDs that are present in the Bray-Curtis distance matrix
in_keep_ids <-  rownames(as.matrix(in_bray_sub))

# Subsetting the metadata so that only those samples are present for analysis in 
# the same order
in_metadata_sub <- gastric_metadata [in_keep_ids, , drop = FALSE]


# Run PERMANOVA Comparing Biological Sex within IN Stage with the subsetted metadata

stats_in <- adonis2(in_bray_sub ~ Gender, data = in_metadata_sub)

stats_in

## Chronic Gastritis (CG) ##

# Select CG gastric disease stage to analyze
cg_stage <- c("Gastric cancer (GC)") 

# Keep ONLY samples from that stage
cg_samples_to_keep = gastric_metadata %>% filter(Group %in% cg_stage) %>% 
  rownames()

# Identify samples within the selected gastric disease stage that are present in 
# the Bray-Curtis distance matrix

cg_dm_ids <- rownames(as.matrix(ps_bray))

cg_samples_to_keep <- gastric_metadata %>%
  filter(Group == cg_stage) %>%
  rownames() %>%
  intersect(cg_dm_ids)

cg_bray_sub <- as.matrix(ps_bray)[cg_samples_to_keep, cg_samples_to_keep] %>% 
  as.dist()

# Obtaining sample IDs that are present in the Bray-Curtis distance matrix
cg_keep_ids <-  rownames(as.matrix(cg_bray_sub))

# Subsetting the metadata so that only those samples are present for analysis in 
# the same order
cg_metadata_sub <- gastric_metadata [cg_keep_ids, , drop = FALSE]


# Run PERMANOVA Comparing Biological Sex within IN Stage with the subsetted metadata

stats_cg <- adonis2(cg_bray_sub ~ Gender, data = cg_metadata_sub)

stats_cg


#### Collection of PERMANOVA Results Across Gastric Disease Stages ####

# Creating table containing p-values obtained from each stage-specific PERMANOVA

permanova_stats <- data.frame(
  Group = c("Healthy control (HC)", "Chronic gastritis (CG)", "Intestinal metaplasia (IM）", 
            "Intraepithelial neoplasia (IN)", "Gastric cancer (GC)"),
  p.value = c(stats_hc$`Pr(>F)`[1],
              stats_cg$`Pr(>F)`[1],
              stats_im$`Pr(>F)`[1],
              stats_in$`Pr(>F)`[1],
              stats_gc$`Pr(>F)`[1])
)


# Extracting p-values from each PERMANOVA result




# Adjust p-values for multiple comparisons using Benjamini-Hochberg (BH) method
permanova_stats <- permanova_stats %>% 
  mutate(padj = p.adjust(p.value, method = "BH"), .after = "p.value")
permanova_stats
