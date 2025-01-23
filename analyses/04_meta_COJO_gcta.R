#!/usr/bin/env Rscript
# Conditional association analysis on meta output

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(data.table)

# Load the dataset
data <- fread("meta_analysis_results.out", sep = "\t", header = TRUE)
filtered_data <- data[n_studies == 2]
head(data)
gcta_input <- filtered_data %>% select(rs_number, reference_allele, other_allele, eaf, beta, se, `p-value`)
gcta_input <- gcta_input %>% rename(SNP = rs_number)
gcta_input <- gcta_input %>% rename(A1 = reference_allele)
gcta_input <- gcta_input %>% rename(A2 = other_allele)
gcta_input <- gcta_input %>% rename(freq = eaf)
gcta_input <- gcta_input %>% rename(b = beta)
gcta_input <- gcta_input %>% rename(se = se)
gcta_input <- gcta_input %>% rename(p = `p-value`)

# Add a new column 'N' with the value 
gcta_input[, N := 622655]

# View the modified data
head(gcta_input)

# Save the updated data to a file

fwrite(gcta_input, "meta_analysis_gwama_two_final1_for_gcta.ma", sep = "\t")

# Running GCTA to identify the lead snps
# Execute the GCTA command
system("./gcta64 --bfile all_phase3 --cojo-file meta_analysis_gwama_two_final1_for_gcta.ma --cojo-slct --out gcta_cojo_results_gwama_two_datasets_all_snps_updated_ma --threads 28")
