# MICB-305-Gastric-Cancer-Research

# Sex Differences in the Gastric Mucosal Microbiome Across Gastric Carcinogenesis

This repository contains the code, analyses, and outputs for a microbiome research project completed for MICB 305: Data Science in Microbiology and Immunology Research (2026W T2) at the University of British Columbia. The project investigated whether biological sex modifies gastric mucosal microbiome composition and predicted function across stages of gastric carcinogenesis.

## Authors

- Christina Xugao
- Ryan Wong
- Parsa Shafiekhani
- Saray Fallahnafari
- Amitis Eskandari

## Project Overview

Gastric cancer develops through a multistage process known as the Correa cascade:

Healthy Control → Chronic Gastritis → Intestinal Metaplasia → Intraepithelial Neoplasia → Gastric Cancer

Recent studies have identified gastric microbiome alterations across disease progression, but the influence of biological sex remains unexplored. This project investigated whether biological sex modifies gastric microbiome composition and predicted function across disease stages.

## Research Question

Does biological sex modify the relationship between gastric disease histopathological stage and gastric microbiome composition and function?

## Objectives

### Aim 1
Determine whether microbial diversity and community structure differ between biological sexes across stages of gastric carcinogenesis.

- Alpha diversity
- Beta diversity

### Aim 2
Determine whether the core gastric mucosal microbiome differs by biological sex.

- Core microbiome analysis

### Aim 3
Determine whether males and females exhibit distinct taxonomic markers across disease stages.

- Indicator taxa analysis
- Differential abundance testing

### Aim 4
Determine whether biological sex influences the predicted functional profile of the gastric mucosal microbiome across disease stages.

- PICRUSt2 functional pathway analysis

## Dataset

This study used a publicly available 16S rRNA sequencing dataset from:

Wang et al. (2020). Changes of the gastric mucosal microbiome associated with histological stages of gastric carcinogenesis. Frontiers in Microbiology. https://doi.org/10.3389/fmicb.2020.00997

### Sample Information

- 310 gastric mucosal biopsy samples
- 132 adult participants
- Northern China cohort

Disease stages:

| Stage | Abbreviation |
|------|-------------|
| Healthy Control | HC |
| Chronic Gastritis | CG |
| Intestinal Metaplasia | IM |
| Intraepithelial Neoplasia | IN |
| Gastric Cancer | GC |

## Methods

- QIIME2 and DADA2 sequence processing
- Taxonomic classification using SILVA database
- Alpha diversity analysis
- Beta diversity analysis
- Core microbiome analysis
- Differential abundance testing (ANCOM-BC2)
- Indicator taxa analysis
- PICRUSt2 functional prediction

## Key Findings

- Sex-associated microbiome differences were stage-specific and most pronounced at the intraepithelial neoplasia stage.
- Females showed higher Faith’s phylogenetic diversity during chronic gastritis.
- Community composition differed by sex only at the intraepithelial neoplasia stage.
- Shared core genera declined across disease progression.
- Brevundimonas and Rhodococcus were identified as male-associated taxa at intraepithelial neoplasia.
- Males showed enrichment of predicted inflammatory and oncogenic pathways at chronic gastritis and intraepithelial neoplasia.

## Repository Structure

### Datasets
Contains the input files used for analysis, including metadata, feature tables, taxonomy assignments, phylogenetic tree files, and processed phyloseq objects.

### R Scripts
Contains scripts used for each major analysis in the project, including:

- Creating Phyloseq object  
- Alpha diversity analysis  
- Beta diversity analysis  
- Core microbiome analysis  
- Differential abundance analysis  
- Indicator taxa analysis  
- Functional pathway analysis  

### Results
Contains outputs generated from the analyses, including:

- **Plots** – final figures  
- **Tables** – indicator genera associated with sex summary  
- **Functional Analysis Annotated pathway** – annotated functional pathway outputs

The plots/ directory contains subdirectories divided by analysis, with generated plots from each stage of the analysis. Many of these are exploratory, intermediate, or unused outputs that were retained for transparency and reproducibility, but are not part of the final manuscript. Formal figures used in the final paper are manually stored separately in Results/Figures/.

### Team Meeting Notes
Notes and records from weekly meetings with our TA, including feedback, discussion points, and project guidance.

### Weekly Meetings 
Planning and scheduling notes for weekly meetings, including progress updates, topics to discuss, completed tasks, and upcoming goals.

