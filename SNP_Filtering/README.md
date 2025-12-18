# SNP Filtering

This folder contains the script used to filter the raw SNP dataset for quality and prepare a high-confidence SNP table for downstream analyses.

## Purpose

We applied several filtering and quality control steps to the raw SNP table to ensure only high-confidence sites were used:

1. Considered sites with sequence coverage ≥20X across all populations.  
2. Removed ~13.8K sites that were not expected to be polymorphic based on founder strain sequences (cf. Phillips et al. 2021).  
3. Removed sites that were not polymorphic in the ancestral population.  
4. Retained sites where the alternate nucleotide frequency fell between 0.02 and 0.98 across the dataset, removing fixed sites and likely sequencing/variant-calling errors.  

These steps produced a high-quality SNP table with **61,281 SNPs**.

## Inputs

- **Raw SNP table:** `SNPtable_raw.txt` (available in the **SNP_Tables** folder on Dryad associated with this project).

## Outputs

- Filtered SNP tables for use in downstream genomic analyses.

