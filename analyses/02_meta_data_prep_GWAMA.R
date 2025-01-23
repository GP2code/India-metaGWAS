#!/usr/bin/env Rscript
# Preps summary statistics for Indian and multi-ancestry datasets, renames columns to standardized meta-analysis fields, and writes out cleaned, tab-separated files. 
# In doing so, it ensures columns like MARKERNAME, NEA, EA, EAF, P, BETA, and SE are properly formatted for GWAMA    

library(tidyverse)
library(dplyr)
library(data.table)

# Preparing the Indian-PD GWAS file for Meta-analysis in GWAMA
# For Indian-PD GWAS
sumstat <- fread("filename")

colnames(sumstat) <- gsub("#", "", colnames(sumstat))

meta_input <- sumstat %>% select(ID, REF, A1, A1_FREQ, P, BETA, SE)
meta_input <- meta_input %>% rename(MARKERNAME = ID)
meta_input <- meta_input %>% rename(NEA = REF)
meta_input <- meta_input %>% rename(EA = A1)
meta_input <- meta_input %>% rename(EAF = A1_FREQ)
meta_input <- meta_input %>% rename(P = P)
meta_input <- meta_input %>% rename(BETA = BETA)
meta_input <- meta_input %>% rename(SE = SE)
fwrite(meta_input, "GWAS1_indian_pd_gwas.txt",sep="\t",row.names=FALSE, quote=FALSE)

#Preparing the Multi-ancestry file after harmonization for Meta-analysis in GWAMA
sumstat <- fread("harmonized_pop_B_multi_ethnic.tsv")

meta_input <- sumstat %>% select(SNP, A2, A1, EAF, P, BETA, SE, CHR, POS)
meta_input <- meta_input %>% rename(MARKERNAME = SNP)
meta_input <- meta_input %>% rename(NEA = A2)
meta_input <- meta_input %>% rename(EA = A1)
meta_input <- meta_input %>% rename(EAF = EAF)
meta_input <- meta_input %>% rename(P = P)
meta_input <- meta_input %>% rename(BETA = BETA)
meta_input <- meta_input %>% rename(SE = SE)

fwrite(meta_input, "GWAS2_multi_ancestry_pd_gwas.txt",sep="\t",row.names=FALSE, quote=FALSE)
