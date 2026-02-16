# In this script, we will run
# 1. Indicator Taxa 
# 2. Differential abundance analysis 

# Load necessary libraries 
library(tidyverse)
library(phyloseq)
library(indicspecies)

# Load object
ps = readRDS('Datasets/phyloseq_taxonomy.rds')


# Indicator species analysis -------------------------

# Aggregate ASVs to the genus level
ps_genus = tax_glom(ps,'Genus')

# convert phyloseq to relative abundance 
ps_relab = transform_sample_counts(ps_genus, function(x) x / sum(x))

#apply abundance filter
ps_filt = filter_taxa(ps_relab, function(x) mean(x) > 0.001, TRUE)
otu_table = data.frame(otu_table(ps_filt))

# Find indicator species for each gastric disease stages - without stratifying
# by biological sex 
set.seed(520)
indval = multipatt(t(otu_table), 
                    cluster= ps_filt@sam_data$Group,
                    control = how(nperm = 999))

summary(indval, indvalcomp = TRUE)

# Extract into table and filter for significance 
indval_table = as.data.frame(indval$sign) %>%
  filter(stat> 0.7 & p.value < 0.05)

# Create dotplot for the indicator taxa separaetd by different groups 
# Get genus names for each taxa ID in ps_filt
genus_map <- tax_table(ps_filt) %>%
  as("matrix") %>%
  as.data.frame() %>%
  rownames_to_column("taxa_id") %>%
  transmute(taxa_id, Genus = as.character(Genus))
# attach taxa_id from rownames + genus
indval_table_genus <- indval_table %>%
  rownames_to_column("taxa_id") %>%
  left_join(genus_map, by = "taxa_id") %>%
  mutate(Genus = if_else(is.na(Genus) | Genus == "", taxa_id, Genus))

# turn dataframe into long format + clean all labels and set stage order 
indval_long <- indval_table_genus  %>%
  pivot_longer(
    cols = starts_with("s."),
    names_to = "Gastric_Disease_Stage",
    values_to = "member"
  ) %>%
  filter(member == 1)%>%
  mutate(
    Genus = str_replace_all(Genus, "g__", ""),
    Gastric_Disease_Stage = substr(Gastric_Disease_Stage, 3, nchar(Gastric_Disease_Stage))) %>%
  mutate(Gastric_Disease_Stage = factor(Gastric_Disease_Stage,
                                        levels = c("Healthy control (HC)", 
                                                   "Chronic gastritis (CG)", 
                                                   "Intestinal metaplasia (IM）", 
                                                   "Intraepithelial neoplasia (IN)", 
                                                   "Gastric cancer (GC)")))

# Keep top N per stage
topN <- 10
plot_df <- indval_long %>%
  group_by(Gastric_Disease_Stage) %>%
  slice_max(stat, n = topN, with_ties = FALSE) %>%
  ungroup()

# Create dotplot 
dotplot_nosex <- ggplot(plot_df, aes(
  x = Gastric_Disease_Stage,
  y = fct_reorder(Genus, stat)
)) +
  geom_point(aes(size = stat, color = Gastric_Disease_Stage,)) +
  scale_size_continuous(range = c(2, 7)) +
  scale_alpha_continuous(range = c(0.4, 1), name = "-log10(p) (capped)") +
  labs(x = NULL, y = "Genus", size = "IndVal stat", color = "Gastric Disease stage") +
  ggtitle("Indicator Genus for Gastric Disease Stage")+
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank()
  )
dotplot_nosex

ggsave("Results/Plots/Ind_Genus_Disease_Stage_nosex.png",
       dotplot_nosex,
       width= 10,
       height = 6)

