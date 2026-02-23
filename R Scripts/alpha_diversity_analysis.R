#Loading Packages
library(tidyverse)
library(phyloseq)


# Loading Object
ps = readRDS("Datasets/phyloseq_taxonomy.rds")

gastric_metadata <- read.delim("Datasets/gastric_cancer_metadata.tsv",
                               row.names = 1)

# Rarefying with a depth of 8000
psrare = ps %>% rarefy_even_depth(sample.size = 8000, rngseed = 421)

#Alpha diversity question: Within each cancer stage, is there a difference in 
#alpha diversity between male and female patients?

#Based on the UJEMI paper by Bains et al. (2025), age was not associated with 
#differences in alpha or beta diversity across any cancer stage. 
#Because of this, I did not control for age in the alpha diversity or 
#beta diversity analyses. However, the paper did find age-related 
#differences in core microbiome composition and predicted functional pathways, 
#especially in young patients. Therefore, age will be controlled for in the 
#core microbiome and functional analyses to avoid confounding.

set.seed(421)
p = plot_richness(psrare,
                  x = "Gender",
                  measures = c("Observed","Shannon","Chao1","Simpson"),
                  color = "Gender")
p

p +
  facet_grid(variable ~ Group, scales = "free_y") +
  theme_bw()

#Extract data
pdata = p$data
str(pdata) # variable = which alpha diversity metric

#value = numeric value of the metric

# Defining a function to plot one alpha diversity metric at a time.
# This avoids repeating the same ggplot code for each metric.

alpha_plot <- function(metric_name){
  
  pdata %>%
    
    # Filtering to only one metric because plotting all metrics together 
    # became messy
    filter(variable == metric_name) %>%
    
    ggplot(aes(x = Gender, y = value, fill = Gender)) +
    
    # Not plotting the default outlier points
    # Plotting all individual points below using jitter
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    
    # Adding jitter to spread points sideways,
    # reducing overlap and making sample distribution clearer.
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    
    # Faceting by cancer stage to compare males vs females within each stage.
    facet_wrap(~Group, nrow = 1) +
    
    # Adding labels and cleaning formatting
    labs(title = metric_name, x = NULL, y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 6.5),
      axis.text.x = element_text(size = 11),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 13, hjust = 0.5)
    )
}

# Wilcoxon Tests for Final Plots of Each Alpha Diversity Metric

# Chao1: Used exact = FALSE because alpha diversity values contain ties.
p_chao1 <- alpha_plot("Chao1") +
  stat_compare_means(
    comparisons = list(c("female","male")),
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    size = 4
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
p_chao1

# Simpson
p_simpson <- alpha_plot("Simpson")+
  stat_compare_means(
    comparisons = list(c("female","male")),
    method = "wilcox.test",
    size = 4
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
p_simpson

# Shannon
p_shannon <- alpha_plot("Shannon") +
  stat_compare_means(
    comparisons = list(c("female","male")),
    method = "wilcox.test",
    size = 4
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
p_shannon

# Observed
p_observed <- alpha_plot("Observed")+
  stat_compare_means(
    comparisons = list(c("female","male")),
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    size = 4
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
p_observed

# Create Alpha_Diversity folder inside Results
dir.create("Results/Plots/Alpha_Diversity")

# Save plots there
ggsave("Results/Plots/Alpha_Diversity/alpha_chao1.png",
       p_chao1, width = 10, height = 4, dpi = 300)

ggsave("Results/Plots/Alpha_Diversity/alpha_simpson.png",
       p_simpson, width = 10, height = 4, dpi = 300)

ggsave("Results/Plots/Alpha_Diversity/alpha_shannon.png",
       p_shannon, width = 10, height = 4, dpi = 300)

ggsave("Results/Plots/Alpha_Diversity/alpha_observed.png",
       p_observed, width = 10, height = 4, dpi = 300)
