#!/bin/bash
# GWAMA to run meta-analysis

# gwama.in; two files

#GWAS1_indian_pd_gwas.txt
#GWAS2_multi_ancestry_pd_gwas.txt

# Fixed effect model
./GWAMA -i gwama.in -o meta_analysis_results -qt

# Random effects model
./GWAMA -i gwama.in --random -o meta_analysis_results -qt