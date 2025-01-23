#!/bin/bash

## Base Data QC
inputdir_basedata='/project/indian-cohort/'
inputdir_targetdata='/project/target-data'
working_dir='/project/..'
output_dir='/project/output'
scripts_dir='/project/scripts/'
prsice_exec='/home/ubuntu/bin/PRSice/'
prefix_bfiles_target='filename_prefix'
basefile='filename'

# BaseData  QC
## Reference for Steps: https://choishingwan.github.io/PRS-Tutorial/base/

## i) Genome build
# - Check that the base and target data are on the same genome build, if not, need to liftover
# base on hg38
gunzip -c $inputdir_basedata$basefile.txt.gz | head 
# Concatenate CHR_BP for SNP ID
gunzip -c $inputdir_basedata$basefile.txt.gz | awk '{OFS="\t"; print $0, $1 "_" $2}' | gzip > $inputdir_basedata$basefile.tsv.gz

## ii) Remove SNPs with imputation INFO > 0.8 and MAF < 0.01
# - SNPs in this base already filtered for MAF > 0.05 and INFO > 0.9
gunzip -c $inputdir_basedata$basefile.tsv.gz | \
awk 'NR==1 || ($11 > 0.01) && ($14 > 0.8) {print}' | gzip  > $inputdir_basedata$basefile.gz

## iii) Remove Duplicate SNPs
gunzip -c $inputdir_basedata$basefile.gz | awk '{seen[$15]++; if(seen[$15]==1){ print}}' | \
gzip  > $inputdir_basedata$basefile.nodup.gz

## iv) Remove Ambiguous SNPs
gunzip -c $inputdir_basedata$basefile.nodup.gz| \
awk '!( ($3=="A" && $4=="T") || \
        ($3=="T" && $4=="A") || \
        ($3=="G" && $4=="C") || \
        ($3=="C" && $4=="G")) {print}' |\
    gzip > $inputdir_basedata$basefile.QC.gz

## v) Format header names
zcat $inputdir_basedata$basefile.QC.gz | awk 'NR==1 {gsub("CHR_BP", "SNP"); gsub("BETA", "OR"); gsub("DR2", "INFO")}1' | \
     gzip > $inputdir_basedata$basefile.base.QC.gz

## keep only required columns
# Select columns
zcat $inputdir_basedata$basefile.base.QC.gz | awk '{print $15, $1, $2, $3, $4, $5, $6, $10}' > $inputdir_basedata$basefile.base.QC.txt 