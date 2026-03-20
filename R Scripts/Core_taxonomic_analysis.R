library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# Load object
ps = readRDS('Datasets/phyloseq_taxonomy.rds')

# Rarefy and convert to relative abundance and use tax_glom to the genus level
psrare = rarefy_even_depth(ps, sample.size = 8000)
ps_rare_relab = transform(psrare, 'compositional')
ps_rare_relab_genus = tax_glom(ps_rare_relab, 'Genus')


# Make core microbiome function ----------------
core_graph <- function(ps, detection, prevalence) {
  # Subset phyloseq object
  sex.male = subset_samples(ps, Gender == 'male')
  sex.female = subset_samples(ps, Gender == 'female')
  
  # Find core members
  ASVs_male = core_members(sex.male, detection=detection, prevalence = prevalence)
  ASVs_female = core_members(sex.female, detection=detection, prevalence = prevalence)
  
  diagram <- ggVennDiagram(
    list(male = ASVs_male, female = ASVs_female),
    set_size = 5
  ) +
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(t =10, r = 100, b = 1, l = 80)
    )
  
  diagram
}

# Run analysis by each stage 

# Subset phyloseq object
sex.male = subset_samples(ps_rare_relab_genus, Gender == 'male')
sex.female = subset_samples(ps_rare_relab_genus, Gender == 'female')

# Find core members
ASVs_male = core_members(sex.male, detection=0.001, prevalence = 0.2)
ASVs_female = core_members(sex.female, detection=0.001, prevalence = 0.2)

diagram <- ggVennDiagram(
  list(male = ASVs_male, female = ASVs_female),
  set_size = 5
) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t =10, r = 100, b = 1, l = 80)
  )

diagram

ggsave("Results/Plots/core_genus_sex.png",
       diagram,
       width = 10,
       height = 8)
