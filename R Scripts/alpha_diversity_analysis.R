# Load libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(broom)
library(writexl)
library(picante)
library(ape)
library(ggsignif)

# Load phyloseq object
ps = readRDS("Datasets/phyloseq_taxonomy.rds")

# Rarefy samples to chosen sequencing depth
psrare = ps %>% rarefy_even_depth(sample.size = 8000, rngseed = 421)
table(sample_sums(psrare))

# Calculate alpha diversity metrics
set.seed(421)
p = plot_richness(
  psrare,
  x = "Gender",
  measures = c("Observed", "Shannon", "Chao1", "Simpson"),
  color = "Gender"
)

# Preview all metrics together
p +
  facet_grid(variable ~ Group, scales = "free_y") +
  theme_bw()

# Extract alpha diversity values for plotting and statistics
pdata = p$data
str(pdata)
# variable = alpha diversity metric
# value = alpha diversity value

# Clean and reorder Group labels
pdata$Group = trimws(as.character(pdata$Group))
pdata$Group[pdata$Group == "Intestinal metaplasia (IM）"] = "Intestinal metaplasia (IM)"
pdata$Group = factor(pdata$Group, levels = group_order)

# Run Wilcoxon tests comparing female vs male within each disease stage
# for each alpha diversity metric, then adjust p-values using BH correction
stats = pdata %>%
  group_by(variable, Group) %>%
  group_modify(~{
    wilcox.test(value ~ Gender, data = .x, conf.int = TRUE) %>%
      broom::tidy()
  }) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

# Find the maximum value for each metric so significance labels can be placed
maxvals = pdata %>%
  group_by(variable) %>%
  summarize(max.y = max(value), .groups = "drop")

# Create a table of significance labels for plotting
annot = stats %>%
  left_join(maxvals, by = "variable") %>%
  mutate(
    ypos = 1.1 * max.y,
    label = ifelse(
      padj > 0.05, "ns",
      ifelse(padj > 0.01, "*",
             ifelse(padj > 0.001, "**", "***"))
    )
  )

# Clean and standardize Group labels in alpha annotation data
annot$Group = trimws(as.character(annot$Group))
annot$Group[annot$Group == "Intestinal metaplasia (IM）"] = "Intestinal metaplasia (IM)"
annot$Group = factor(annot$Group, levels = group_order)

# Function to plot one alpha diversity metric at a time
alpha_plot <- function(metric_name) {
  
  plot_data <- pdata %>%
    filter(variable == metric_name)
  
  annot_data <- annot %>%
    filter(variable == metric_name)
  
  ggplot(plot_data, aes(x = Gender, y = value, fill = Gender)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1.5) +
    facet_wrap(~Group, nrow = 1) +
    geom_signif(
      data = annot_data %>% filter(label %in% c("*", "**", "***")),
      aes(
        xmin = "female",
        xmax = "male",
        annotations = label,
        y_position = ypos),
      manual = TRUE,
      inherit.aes = FALSE,
      textsize = 6,
      tip_length = 0.02,
      vjust = 0.3) +
    geom_signif(
      data = annot_data %>% filter(label == "ns"),
      aes(xmin = "female",
        xmax = "male",
        annotations = label,
        y_position = ypos),
      manual = TRUE,
      inherit.aes = FALSE,
      textsize = 4,
      tip_length = 0.02,
      vjust = 0.3
    ) +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(labels = c("female" = "Female", "male" = "Male")) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 6.5),
      axis.text.x = element_text(size = 8, angle = 50, hjust = 1),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 13, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.115)))
}

# Generate one plot for each alpha diversity metric
p_observed = alpha_plot("Observed") + labs(y = "Observed ASVs")
p_shannon  = alpha_plot("Shannon")  + labs(y = "Shannon Diversity Index")
p_chao1    = alpha_plot("Chao1")    + labs(y = "Chao1 Richness")
p_simpson  = alpha_plot("Simpson")  + labs(y = "Simpson Diversity Index")


# Save alpha diversity plots
ggsave("Results/Plots/Alpha_Observed.png",
       plot = p_observed, height = 4, width = 12)

ggsave("Results/Plots/Alpha_Shannon.png",
       plot = p_shannon, height = 4, width = 12)

ggsave("Results/Plots/Alpha_Chao1.png",
       plot = p_chao1, height = 4, width = 12)

ggsave("Results/Plots/Alpha_Simpson.png",
       plot = p_simpson, height = 4, width = 12)

# Save alpha diversity statistics
write_xlsx(
  list("Alpha Diversity Stats" = stats),
  "../Results/Tables/Alpha_Diversity_Stats.xlsx"
)


# Faith's Phylogenetic Diversity


# Convert OTU table to a matrix and make sure taxa are columns
otu = as(otu_table(psrare), "matrix")
if (taxa_are_rows(psrare) == TRUE) otu <- t(otu)

# Extract the phylogenetic tree from the phyloseq object
tree = phy_tree(psrare)

# Make sure OTU table taxa names match tree tip labels
stopifnot(all(colnames(otu) %in% tree$tip.label))
tree = drop.tip(tree, setdiff(tree$tip.label, colnames(otu)))
tree = keep.tip(tree, colnames(otu))

# Root the tree if needed
if (!is.rooted(tree)) {
  tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
}

# Convert counts to presence/absence for Faith's PD
otu_pa = otu
otu_pa[otu_pa > 0] = 1

# Calculate Faith's Phylogenetic Diversity
faith = pd(otu_pa, tree, include.root = TRUE) %>%
  select(-SR)

# Add sample IDs and metadata for plotting
faith_df = faith %>%
  rownames_to_column("SampleID")

meta_df = data.frame(sample_data(psrare)) %>%
  rownames_to_column("SampleID")

faith_df = left_join(faith_df, meta_df, by = "SampleID")

# Clean and reorder Group labels
faith_df$Group = trimws(as.character(faith_df$Group))
faith_df$Group[faith_df$Group == "Intestinal metaplasia (IM）"] = "Intestinal metaplasia (IM)"

# Order groups by gastric cancer progression
group_order = c("Healthy control (HC)",
                "Chronic gastritis (CG)",
                "Intestinal metaplasia (IM)",
                "Intraepithelial neoplasia (IN)",
                "Gastric cancer (GC)")

faith_df$Group = factor(faith_df$Group, levels = group_order)

# Run Wilcoxon tests comparing female vs male within each disease stage
faith_stats = faith_df %>%
  group_by(Group) %>%
  group_modify(~{
    wilcox.test(PD ~ Gender, data = .x, conf.int = TRUE) %>%
      broom::tidy()
  }) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH"))

# Use one common height for significance labels across the full Faith's PD plot
faith_max = max(faith_df$PD, na.rm = TRUE)

faith_annot = faith_stats %>%
  mutate(
    ypos = 1.08 * faith_max,
    label = ifelse(
      padj > 0.05, "ns",
      ifelse(padj > 0.01, "*",
             ifelse(padj > 0.001, "**", "***"))
    ),
    label_size = ifelse(label == "*", 8, 6)
  )

# Clean and standardize Group labels in faith annotation data
faith_annot$Group = trimws(as.character(faith_annot$Group))
faith_annot$Group[faith_annot$Group == "Intestinal metaplasia (IM）"] = "Intestinal metaplasia (IM)"
faith_annot$Group = factor(faith_annot$Group, levels = group_order)

# Create the Faith's PD plot
p_faith = ggplot(faith_df, aes(x = Gender, y = PD, fill = Gender)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.5) +
  facet_wrap(~Group, nrow = 1) +
  geom_signif(data = faith_annot %>% filter(label == "*"),
              aes(xmin = "female", xmax = "male",
                  annotations = label,
                  y_position = ypos),
              manual = TRUE,
              inherit.aes = FALSE,
              textsize = 6,   
              tip_length = 0.02,
              vjust = 0.3) +
  geom_signif(data = faith_annot %>% filter(label == "ns"),
              aes(xmin = "female", xmax = "male",
                  annotations = label,
                  y_position = ypos),
              manual = TRUE,
              inherit.aes = FALSE,
              textsize = 4,   
              tip_length = 0.02,
              vjust = 0.3) +
  scale_size_identity() +
  labs(y = "Phylogenetic Diversity (Faith's PD)", x = NULL) +
  scale_x_discrete(labels = c("female" = "Female", "male" = "Male")) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 6.5),
    axis.text.x = element_text(size = 8, angle = 50, hjust = 1),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 13, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0.02, 0.115)))

p_faith

# Save Faith's PD plot
ggsave("Results/Plots/Alpha_Faith_PD.png",
       plot = p_faith, height = 4, width = 12)

# Save Faith's PD statistics
write_xlsx(
  list("Faith PD Stats" = faith_stats),
  "../Results/Tables/Faith_PD_Stats.xlsx"
)
