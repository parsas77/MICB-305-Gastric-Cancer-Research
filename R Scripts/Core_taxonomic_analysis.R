library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# Load object
ps = readRDS('Datasets/phyloseq_taxonomy.rds')

# Rarefy if desired and convert to relative abundance, use tax_glom (optional)
psrare = rarefy_even_depth(ps, sample.size = 8000)
ps_rare_relab = transform(psrare, 'compositional')
ps_rare_relab_genus = tax_glom(ps_rare_relab, 'Genus')

# Subset phyloseq object
sex.male = subset_samples(ps_rare_relab_genus, Gender == 'male')
sex.female = subset_samples(ps_rare_relab_genus, Gender == 'female')

# Find core members
ASVs_male = core_members(sex.male, detection=0.001, prevalence = 0.8)
ASVs_female = core_members(sex.female, detection=0.001, prevalence = 0.8)

diagram<- ggVennDiagram(list(ASVs_male, ASVs_female),
              set_size = 5,
              category.names = c('male','female'))

ggsave("Results/Plots/core_genus_sex.png",
       diagram,
       width = 8,
       height = 5)
# 1 core taxa unique to male and no unique core taxa for female so I 
# don't think we need to go futher into each stages. 