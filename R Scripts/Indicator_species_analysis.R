# In this script, we will run Indicator Taxa Analysis

# Load necessary libraries 
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(gt)

# Load object
ps = readRDS('Datasets/phyloseq_taxonomy.rds')

# Aggregate ASVs to the genus level
# convert phyloseq to relative abundance 
#apply abundance filter
ps_genus = tax_glom(ps,'Genus')
ps_relab = transform_sample_counts(ps_genus, function(x) x / sum(x))
ps_filt = filter_taxa(ps_relab, function(x) mean(x) > 0.001, TRUE)
otu_table = data.frame(otu_table(ps_filt))

# Indicator species analysis --------------
# Stratified by Gastric Disease Stage only -------------------------

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
       height = 10)

# Individual analysis for each sex between stages  -------------------------

# Function 1 ------------------

run_indval <- function(ps_obj, otu_table, seed = 520, cat) {
  set.seed(520)
  indval = multipatt(t(otu_table), 
                     cluster= ps_obj@sam_data$Group,
                     control = how(nperm = 999))
  
  indval_table = as.data.frame(indval$sign) %>%
    filter(stat> 0.7 & p.value < 0.05)
  genus_map <- tax_table(ps_obj) %>%
    as("matrix") %>%
    as.data.frame() %>%
    rownames_to_column("taxa_id") %>%
    transmute(taxa_id, Genus = as.character(Genus))
  # attach taxa_id from rownames + genus
  indval_table_genus <- indval_table %>%
    rownames_to_column("taxa_id") %>%
    left_join(genus_map, by = "taxa_id") %>%
    mutate(Genus = if_else(is.na(Genus) | Genus == "", taxa_id, Genus))
  
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

  plot_df <- indval_long %>%
    group_by(Gastric_Disease_Stage) %>%
    slice_max(stat, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  dotplot <- ggplot(plot_df, aes(
    x = Gastric_Disease_Stage,
    y = fct_reorder(Genus, stat)
  )) +
    geom_point(aes(size = stat, color = Gastric_Disease_Stage,)) +
    scale_size_continuous(range = c(2, 7)) +
    scale_alpha_continuous(range = c(0.4, 1), name = "-log10(p) (capped)") +
    labs(x = NULL, y = "Genus", size = "IndVal stat", color = "Gastric Disease stage") +
    ggtitle(paste("Indicator Genus for Gastric Disease Stage", "-", cat))+
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.minor = element_blank()
    )
  dotplot
}


ps_f <- subset_samples(ps_filt, Gender == "female")
otu_table_f = data.frame(otu_table(ps_f))

dotplot <- run_indval(ps_f, 
                      otu_table_f, 
                      seed = 520, 
                      cat = "Female") 
dotplot
ggsave("Results/Plots/Ind_Genus_Disease_Stage_Female.png",
       dotplot,
       width= 10,
       height = 10)

ps_m <- subset_samples(ps_filt, Gender == "male")
otu_table_m = data.frame(otu_table(ps_m))
dotplot <- run_indval(ps_f, 
                      otu_table_f, 
                      seed = 520, 
                      cat = "Male") 
dotplot
ggsave("Results/Plots/Ind_Genus_Disease_Stage_Male.png",
       dotplot,
       width= 10,
       height = 10)


# Comparison between sex within each stage  -------------------------
# No significant indicator taxa between male and female for all 4 stages except for IN

# Function 2 ---------------------
run_indval <- function(ps_obj, otu_table, seed = 520, cat) {
  set.seed(520)
  indval = multipatt(t(otu_table), 
                     cluster= ps_obj@sam_data$Gender,
                     control = how(nperm = 999))
  
  indval_table = as.data.frame(indval$sign) %>%
    filter(stat> 0.7 & p.value < 0.05)
  genus_map <- tax_table(ps_obj) %>%
    as("matrix") %>%
    as.data.frame() %>%
    rownames_to_column("taxa_id") %>%
    transmute(taxa_id, Genus = as.character(Genus))
  # attach taxa_id from rownames + genus
  indval_table_genus <- indval_table %>%
    rownames_to_column("taxa_id") %>%
    left_join(genus_map, by = "taxa_id") %>%
    mutate(Genus = if_else(is.na(Genus) | Genus == "", taxa_id, Genus))
  
  indval_long <- indval_table_genus  %>%
    pivot_longer(
      cols = starts_with("s."),
      names_to = "Biological_Sex",
      values_to = "member"
    ) %>%
    filter(member == 1)%>%
    mutate(
      Genus = str_replace_all(Genus, "g__", ""),
      Biological_Sex = substr(Biological_Sex, 3, nchar(Biological_Sex))) %>%
    mutate(Biological_Sex = factor(Biological_Sex,
                                   levels = c("male",
                                              "female")))
  
  plot_df <- indval_long %>%
    group_by(Biological_Sex) %>%
    slice_max(stat, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  dotplot <- ggplot(plot_df, aes(
    x = Biological_Sex,
    y = fct_reorder(Genus, stat)
  )) +
    geom_point(aes(size = stat, color = Biological_Sex,)) +
    scale_size_continuous(range = c(2, 7)) +
    scale_alpha_continuous(range = c(0.4, 1), name = "-log10(p) (capped)") +
    labs(x = NULL, y = "Genus", size = "IndVal stat", color = "Biological Sex") +
    ggtitle(paste("Indicator Genus for Sex", "-", cat))+
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.minor = element_blank()
    )
  dotplot
}


hc <- subset_samples(ps_filt, Group == "Healthy control (HC)")
otu_table_hc = data.frame(otu_table(hc))

dotplot <- run_indval(hc, 
                      otu_table_hc, 
                      seed = 520, 
                      cat = "Healthy control") 
dotplot

CG <- subset_samples(ps_filt, Group == "Chronic gastritis (CG)")
otu_table_CG = data.frame(otu_table(CG))

dotplot <- run_indval(CG, 
                      otu_table_CG, 
                      seed = 520, 
                      cat = "Chronic gastritis (CG)") 
dotplot

IM <- subset_samples(ps_filt, Group == "Intestinal metaplasia (IM）")
otu_table_IM = data.frame(otu_table(IM))

dotplot <- run_indval(IM, 
                      otu_table_IM, 
                      seed = 520, 
                      cat = "Intestinal metaplasia") 
dotplot

IN <- subset_samples(ps_filt, Group == "Intraepithelial neoplasia (IN)")
otu_table_IN = data.frame(otu_table(IN))

dotplot <- run_indval(IN, 
                      otu_table_IN, 
                      seed = 520, 
                      cat = "Intraepithelial neoplasia (IN)") 
dotplot

ggsave("Results/Plots/IN_Ind_Genus_Disease_Stage.png",
       dotplot,
       width= 10,
       height = 10)

GC <- subset_samples(ps_filt, Group == "Gastric cancer (GC)")
otu_table_GC = data.frame(otu_table(GC))

dotplot <- run_indval(GC, 
                      otu_table_GC, 
                      seed = 520, 
                      cat = "Gastric cancer (GC)") 
dotplot

# Summary table for male at IN stage ----------------
# Since there are two indicator taxa unqiue to male at the IN stage, I will
# make a table displaying the two genus and the associated stats 

# Re-run the analysis for IN stage 
set.seed(520)
indval = multipatt(t(otu_table_IN), 
                   cluster= IN@sam_data$Gender,
                   control = how(nperm = 999))

indval_table = as.data.frame(indval$sign) %>%
 filter(stat> 0.7 & p.value < 0.05)

genus_map <- tax_table(IN) %>%
  as("matrix") %>%
  as.data.frame() %>%
  rownames_to_column("taxa_id") %>%
  transmute(taxa_id, Genus = as.character(Genus))

indval_table_genus <- indval_table %>%
  rownames_to_column("taxa_id") %>%
  left_join(genus_map, by = "taxa_id") %>%
  mutate(
    Genus = if_else(is.na(Genus) | Genus == "", taxa_id, Genus),
    Genus = str_remove(Genus, "^g__"),
    Genus = if_else(Genus == "" | is.na(Genus), "Unclassified", Genus),
    Biological_Sex = case_when(
      s.female == 1 & s.male == 0 ~ "Female",
      s.female == 0 & s.male == 1 ~ "Male",
      s.female == 1 & s.male == 1 ~ "Female and male",
      TRUE ~ "Other"
    ),
    Indicator_value = round(stat, 3),
    p_value = if_else(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  ) %>%
  select(Genus, Biological_Sex, Indicator_value, p_value) %>%
  arrange(Biological_Sex, desc(Indicator_value))
indval_table_genus

indicator_pub_table <- indval_table_genus %>%
  gt() %>%
  tab_header(
    title = md("**Indicator genera associated with sex**"),
    subtitle = "Genera with indicator value > 0.7 and p < 0.05"
  ) %>%
  cols_label(
    Genus = "Genus",
    Biological_Sex = "Biological Sex",
    Indicator_value = "Indicator value",
    p_value = md("**p**-value")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = Genus)
  )%>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 13,
    heading.subtitle.font.size = 11,
    data_row.padding = px(5),
    table.border.top.width = px(1.5),
    table.border.bottom.width = px(1.5),
    column_labels.border.bottom.width = px(1),
    row.striping.background_color = "grey95"
  ) %>%
  opt_row_striping() %>%
  tab_source_note(
    source_note = "Indicator taxa were identified using multipatt analysis with 999 permutations."
  )

indicator_pub_table

# Save the table
gtsave(indicator_pub_table, "Results/Tables/IN_indicator_genera_sex.png")
