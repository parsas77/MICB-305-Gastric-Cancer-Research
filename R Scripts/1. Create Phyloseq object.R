library(tidyverse)
library(phyloseq)

# Load datasets
taxonomy = read.delim('Datasets/taxonomy.tsv', row.names = 1)
tree = read_tree('Datasets/tree.nwk')
counts = read.delim('Datasets/feature-table.txt', skip=1, row.names=1) 
metadata = read.delim('Datasets/gastric_cancer_metadata.tsv', row.names = 1) 

# Wrangle Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Taxonomy
taxonomy_formatted = taxonomy %>% 
  separate(col = Taxon, 
           into = c('Domain','Phylum','Class','Order',
                    'Family','Genus','Species'),
           sep=';', fill='right') %>% 
  select(-Confidence) %>% 
  as.matrix()

# Counts
counts_formatted = counts %>% as.matrix()

# Metadata all cleaned and does not require further wrangling 

# Create the phyloseq object
ps = phyloseq(sample_data(metadata),
              otu_table(counts_formatted, taxa_are_rows = T),
              tax_table(taxonomy_formatted),
              tree)

# Save as .rds or .Rdata object  
saveRDS(ps,'Datasets/phyloseq_taxonomy.rds')
