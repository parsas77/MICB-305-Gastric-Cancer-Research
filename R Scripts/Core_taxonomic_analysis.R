library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(patchwork)

# Load object
ps = readRDS('Datasets/phyloseq_taxonomy.rds')

# Rarefy and convert to relative abundance and use tax_glom to the genus level
psrare = rarefy_even_depth(ps, sample.size = 8000)
ps_rare_relab = transform(psrare, 'compositional')
ps_rare_relab_genus = tax_glom(ps_rare_relab, 'Genus')

# Make core microbiome function ----------------

get_stage_fill_max <- function(ps, detection = 0.001, prevalence = 0.2) {
  sex.male <- subset_samples(ps, Gender == "male")
  sex.female <- subset_samples(ps, Gender == "female")
  
  ASVs_male <- core_members(sex.male, detection = detection, prevalence = prevalence)
  ASVs_female <- core_members(sex.female, detection = detection, prevalence = prevalence)
  
  venn_obj <- Venn(list(Female = ASVs_female, Male = ASVs_male))
  venn_dat <- process_data(venn_obj)
  
  max(venn_region(venn_dat)$count, na.rm = TRUE)
}

# Find maximum value for each stage to set limit for scale fill 
hc <- subset_samples(ps_rare_relab_genus, Group == "Healthy control (HC)")
CG <- subset_samples(ps_rare_relab_genus, Group == "Chronic gastritis (CG)")
IM <- subset_samples(ps_rare_relab_genus, Group == "Intestinal metaplasia (IM）")
IN <- subset_samples(ps_rare_relab_genus, Group == "Intraepithelial neoplasia (IN)")
GC <- subset_samples(ps_rare_relab_genus, Group == "Gastric cancer (GC)")

scale_limit <- max(get_stage_fill_max(hc),
                   get_stage_fill_max(CG), 
                   get_stage_fill_max(IM), 
                   get_stage_fill_max(IN), 
                   get_stage_fill_max(GC))


core_graph <- function(ps, detection, prevalence, title = NULL) {
  # Subset phyloseq object
  sex.male <- subset_samples(ps, Gender == "male")
  sex.female <- subset_samples(ps, Gender == "female")
  
  # Find core members
  ASVs_male <- core_members(sex.male, detection = detection, prevalence = prevalence)
  ASVs_female <- core_members(sex.female, detection = detection, prevalence = prevalence)
  
  ggVennDiagram(
    list(Female = ASVs_female, Male = ASVs_male),
    set_size = 8,
    label_size = 12
  ) +
    coord_cartesian(clip = "off") +
    ggtitle(title) +
    scale_fill_continuous(limits = c(0, scale_limit))+
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      plot.margin = margin(t =10, r = 10, b = 1, l = 10),
      panel.border = element_blank()
    )
}

# Run analysis by each stage -------------- 
hc_plot <- core_graph(hc, 0.001, 0.2, "Healthy control")+
  theme(legend.position = "none")
hc_plot 

ggsave("Results/Plots/core_microbiome/core_genus_sex_HC.png",
       hc_plot,
       width = 13,
       height = 13)

CG_plot <- core_graph(CG, 0.001, 0.2, "Chronic gastritis")+
  theme(legend.position = "none")
CG_plot
ggsave("Results/Plots/core_microbiome/core_genus_sex_CG.png",
       CG_plot,
       width = 13,
       height = 13)

IM <- subset_samples(ps_rare_relab_genus, Group == "Intestinal metaplasia (IM）")
IM_plot <- core_graph(IM, 0.001, 0.2, "Intestinal metaplasia")
IM_plot
ggsave("Results/Plots/core_microbiome/core_genus_sex_IM.png",
       IM_plot,
       width = 13,
       height = 13)


IN_plot <- core_graph(IN, 0.001, 0.2, "Intraepithelial neoplasia")+
  theme(legend.position = "none")
IN_plot
ggsave("Results/Plots/core_microbiome/core_genus_sex_IN.png",
       IN_plot,
       width = 13,
       height = 13)


GC_plot <- core_graph(GC, 0.001, 0.2, "Gastric cancer")+
  theme(legend.position = "none")
GC_plot
ggsave("Results/Plots/core_microbiome/core_genus_sex_GC.png",
       GC_plot,
       width = 13,
       height = 13)

# Make combined graphs

combined_plot <- wrap_plots(
  hc_plot, CG_plot, IM_plot,
  IN_plot, GC_plot,
  ncol = 3
)
combined_plot

ggsave("Results/Plots/core_microbiome/core_genus_combined.png",
       combined_plot,
       width = 26,
       height = 18)


