# Load libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(broom)
library(writexl)
library(picante)
library(ape)

# Load phyloseq object
ps = readRDS("../Datasets/phyloseq_taxonomy.rds")

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
      padj > 0.05, "NS",
      ifelse(padj > 0.01, "*",
             ifelse(padj > 0.001, "**", "***"))
    )
  )

# Function to plot one alpha diversity metric at a time
alpha_plot <- function(metric_name) {
  
  plot_data <- pdata %>%
    filter(variable == metric_name)
  
  annot_data <- annot %>%
    filter(variable == metric_name)
  
  ggplot(plot_data, aes(x = Gender, y = value, fill = Gender)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    facet_wrap(~Group, nrow = 1) +
    geom_text(
      data = annot_data,
      aes(x = 1.5, y = ypos, label = label),
      inherit.aes = FALSE,
      size = 5
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 6.5),
      axis.text.x = element_text(size = 11),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 13, hjust = 0.5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))
}

# Generate one plot for each alpha diversity metric
p_observed = alpha_plot("Observed") + labs(y = "Observed ASVs")
p_shannon  = alpha_plot("Shannon")  + labs(y = "Shannon Diversity Index")
p_chao1    = alpha_plot("Chao1")    + labs(y = "Chao1 Richness")
p_simpson  = alpha_plot("Simpson")  + labs(y = "Simpson Diversity Index")

# Save alpha diversity plots
ggsave("../Results/Plots/Alpha_Observed.jpeg",
       plot = p_observed, height = 4, width = 12)

ggsave("../Results/Plots/Alpha_Shannon.jpeg",
       plot = p_shannon, height = 4, width = 12)

ggsave("../Results/Plots/Alpha_Chao1.jpeg",
       plot = p_chao1, height = 4, width = 12)

ggsave("../Results/Plots/Alpha_Simpson.jpeg",
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
      padj > 0.05, "NS",
      ifelse(padj > 0.01, "*",
             ifelse(padj > 0.001, "**", "***"))
    ),
    label_size = ifelse(label == "NS", 6, 8)
  )

# Create the Faith's PD plot
p_faith = ggplot(faith_df, aes(x = Gender, y = PD, fill = Gender)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.5) +
  facet_wrap(~Group, nrow = 1) +
  geom_text(
    data = faith_annot,
    aes(x = 1.5, y = ypos, label = label, size = label_size),
    inherit.aes = FALSE
  ) +
  scale_size_identity() +
  labs(y = "Phylogenetic Diversity (Faith's PD)", x = NULL) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

p_faith

# Save Faith's PD plot
ggsave("../Results/Plots/Alpha_Faith_PD.jpeg",
       plot = p_faith, height = 4, width = 12)

# Save Faith's PD statistics
write_xlsx(
  list("Faith PD Stats" = faith_stats),
  "../Results/Tables/Faith_PD_Stats.xlsx"
)
