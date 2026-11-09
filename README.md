# Medina River Microbial Community and Antibiotic Resistance Gene Analysis
Overview

This repository contains Python and R scripts used for the downstream analysis and visualization of microbial community composition and antimicrobial resistance genes (ARGs) in water samples collected from the Medina River watershed, Texas, USA.
The scripts support the analyses presented in a manuscript submitted to Environmental Microbiome.
The repository includes R and python scripts for microbial community diversity and taxonomy analyses, and ARG analyses, including ARG diversity and prevalence, antibiotic resistance profiles, ARG-host associations for river water samples.

Repository Contents:

Metadata:
Metadata.csv - Sample metadata used to associate analytical results with sampling locations and collection events.

Raw sequencing data:

Raw sequencing data generated for this study are publicly available through the NCBI Sequence Read Archive (SRA).

1. 16S rRNA gene amplicon sequencing data: NCBI BioProject accession: PRJNA1442845
2. Nanopore metagenomic sequencing data: NCBI BioProject accession: PRJNA1454989

Python scripts:
1. argo_amr_summary.py - Summarizes ARG richness, diversity, prevalence, drug-class profiles, plasmid-associated ARGs, and ARG-carrying taxa from ARGO outputs.
2. argo_lineages_cpg_subtype_samples.py - Processes ARG lineage and subtype information and calculates or summarizes ARG abundance normalized as copies per Gbp.
3. generate_ARG_host_comparison_tables.py - Generates comparison tables for ARG-host associations.
4. summarize_ARG_carrier_host_associations.py - Summarizes associations between ARGs and their predicted microbial hosts.
5. summarize_ARG_drug_class.py - Summarizes ARGs by antibiotic resistance class.

R scripts:
1. R_script-bacterial_community_beta-diversity.txt - Analyzes and visualizes differences in bacterial community composition among samples.
2. R_script-ARG_braycurtis_distance_heatmap.txt - Generates a heatmap of Bray–Curtis dissimilarities based on ARG profiles.
3. R_script-ARG_braycurtis_distance_NMDS.txt - Performs non-metric multidimensional scaling (NMDS) analysis of ARG profiles using Bray–Curtis dissimilarities.
4. R_script-ARG_prevalence_plot.txt - Visualizes ARG prevalence across samples.
5. R_script-Clinic_ARG_profile_heatmap.txt - Generates a heatmap of ARG profiles associated with clinically relevant resistance categories.
6. R_script-Drug_class_profile_heatmap.txt - Generates a heatmap showing ARG abundance across antibiotic resistance classes.
7. R_script-Temporal_ARGs_abundance_diversity_plasmid.txt - Visualizes temporal patterns in ARG abundance, diversity, and plasmid association.

