# Load necessary libraries 
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

# Load object
ps = readRDS('Datasets/phyloseq_taxonomy.rds')

# Differential abundance analysis -------------------------

# Aggregate ASVs to the genus level
ps_genus = tax_glom(ps,'Genus')


# Fix a formatting issue in Group label (IM had a different parenthesis character)
ps_genus@sam_data$Group <- case_when(
  ps_genus@sam_data$Group == "Intestinal metaplasia (IM）" ~ "Intestinal metaplasia (IM)",
  TRUE ~ ps_genus@sam_data$Group
)


# Make Gender a factor with female as the reference
ps_genus@sam_data$Gender = factor(ps_genus@sam_data$Gender,
                                 levels = c('female',
                                            'male'))

# Make Group a factor
ps_genus@sam_data$Group = factor(ps_genus@sam_data$Group,
                                 levels = c("Healthy control (HC)", 
                                            "Chronic gastritis (CG)", 
                                            "Intestinal metaplasia (IM)", 
                                            "Intraepithelial neoplasia (IN)", 
                                            "Gastric cancer (GC)"))


# 1) Healthy only: Male vs Female
ps_healthy <- subset_samples(ps_genus, Group == "Healthy control (HC)")

set.seed(522)
healthy_out <- ancombc2(
  data = ps_healthy,
  fix_formula = "Gender",
  p_adj_method = "BH",
  prv_cut = 0.1
)

healthy_statistical_table = healthy_out$res

## NO SIGNIFANCE BETWEEN MALE AND FEMALE HEALTHY GROUPS



# 2) Chronic gastritis only: Male vs Female
ps_chronic_gastritis <- subset_samples(ps_genus, Group == "Chronic gastritis (CG)")

set.seed(522)
chronic_gastritis_out <- ancombc2(
  data = ps_chronic_gastritis,
  fix_formula = "Gender",
  p_adj_method = "BH",
  prv_cut = 0.1
)

chronic_gastritis_statistical_table = chronic_gastritis_out$res

## NO SIGNIFANCE BETWEEN MALE AND FEMALE CG GROUPS



# 3) Intestinal metaplasia only: Male vs Female
ps_intestinal_metaplasia <- subset_samples(ps_genus, Group == "Intestinal metaplasia (IM)")

set.seed(522)
intestinal_metaplasia_out <- ancombc2(
  data = ps_intestinal_metaplasia,
  fix_formula = "Gender",
  p_adj_method = "BH",
  prv_cut = 0.1
)

intestinal_metaplasia_statistical_table = intestinal_metaplasia_out$res

## NO SIGNIFANCE BETWEEN MALE AND FEMALE IM GROUPS



# 4) Intraepithelial neoplasia only: Male vs Female
ps_intraepithelial_neoplasia <- subset_samples(ps_genus, Group == "Intraepithelial neoplasia (IN)")

set.seed(522)
intraepithelial_neoplasia_out <- ancombc2(
  data = ps_intraepithelial_neoplasia,
  fix_formula = "Gender",
  p_adj_method = "BH",
  prv_cut = 0.1
)

intraepithelial_neoplasia_statistical_table = intraepithelial_neoplasia_out$res

# Filter stats to only include taxa that are differentially abundant between male and female within (IN)
intraepithelial_neoplasia_taxa <- intraepithelial_neoplasia_statistical_table %>%
  filter(diff_robust_Gendermale == TRUE)

# Plot the log 2 fold changes of differential abundance of intraepithelial neoplasia (male relative
# to female) 
intraepithelial_neoplasia_plot <- intraepithelial_neoplasia_taxa%>%
  mutate(
    Genus = as.character(tax_table(ps_genus)[taxon, "Genus"]),
    Genus = str_replace_all(Genus, "g__", "")
  ) %>%
  ggplot(aes(Genus, lfc_Gendermale, fill = Genus)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Bacterial Genus",
    y = "Log Fold Change (Male relative to Female)",
    title = "Differential Abundance (Intraepithelial Neoplasia)"
  )

intraepithelial_neoplasia_plot


# 5) Gastric cancer only: Male vs Female
ps_gastric_cancer <- subset_samples(ps_genus, Group == "Gastric cancer (GC)")

set.seed(522)
gastric_cancer_out <- ancombc2(
  data = ps_gastric_cancer,
  fix_formula = "Gender",
  p_adj_method = "BH",
  prv_cut = 0.1
)

gastric_cancer_statistical_table = gastric_cancer_out$res

## NO SIGNIFANCE BETWEEN MALE AND FEMALE GC GROUPS


# Save plot to GitHub results folder
ggsave("Results/Plots/IN_Genus_Differential_Abundance_Male_vs_Female.png",
       intraepithelial_neoplasia_plot,
       width = 10,
       height = 6)


