#!/bin/bash

## Target Data QC
# - Imputed variants for Luxgiant data
SCRIPTS_DIR='/project/scripts'
inputdir_basedata='/project/indian-cohort/'
inputdir_targetdata='/project/target-data/'
working_dir='/project/..'
output_dir='/project/output'
prsice_exec='/home/ubuntu/bin/PRSice/'
prefix_bfiles_target='filename_prefix'
basefile='filename'
LDregions_file=/project/high-LD-regions.txt 
threadnum=50
memory_num=300000


# A) TargetData  QC
## Reference for Steps: https://choishingwan.github.io/PRS-Tutorial/target/

## Genome build
# - Check that the base and target data are on the same genome buil, if not, need to liftover
# Target on hg38

## Standard GWAS QC
# Collect  failed SNPs in a list
plink --bfile $inputdir_targetdata$prefix_bfiles_target \
      --chr 1-22 \
      --maf 0.01 \
      --geno 0.01 \
      --hwe 1e-6 \
      --mind 0.01 \
      --memory ${memory_num} \
      --threads ${threadnum} \
      --make-bed \
      --out $output_dir$prefix_bfiles_target-target-qced

## Substitute rs_IDs with CHR_POS where there are no
mv $output_dir$prefix_bfiles_target-target-qced.bim $woutput_dir$prefix_bfiles_target-qced.bim.old
awk '{OFS="\t"; print $1, $1 "_" $4, $3, $4, $5, $6}' $output_dir$prefix_bfiles_target-qced.bim.old > $output_dir$prefix_bfiles_target-target-qced.bim 


# B) LD Prune data
## a. Excluding high-LD and HLA regions
plink --bfile $output_dir$prefix_bfiles_target-target-qced \
      --exclude $LDregions_file \
      --indep-pairwise 50 5 0.2 \
      --memory 200000 \
      --threads ${threadnum} \
      --make-bed \
      --out $output_dir$prefix_bfiles_target-pruning \

## b. LD prune (LD -> highly non random correlated blocks of SNPs)
plink --bfile $output_dir$prefix_bfiles_target-pruning \
      --extract $output_dir$prefix_bfiles_target-pruning.prune.in \
      --make-bed \
      --out $output_dir$prefix_bfiles_target-LDpruned \
      --memory ${memory_num} \
      --threads ${threadnum}

## c. Generate PCA with prunned data
plink --bfile $output_dirr$prefix_bfiles_target-LDpruned \
      --pca 10 \
      --memory ${memory_num} \
      --threads ${threadnum} \
      --out $output_dir$prefix_bfiles_target-grm-pca
      
## d.Phenotype 
# Missing coded as NA or -9 for binary traits and NA for quantitative traits
#Extract phenotype from .fam file
awk '{print $1, $2, $6 }' $output_dir$prefix_bfiles_target-target-qced.fam > $output_dir$prefix_bfiles_target-pheno.txt
# Add column names to phentype file
awk 'BEGIN {print "FID IID PD"} 1' $output_dir$prefix_bfiles_target-pheno.txt > $output_dir$prefix_bfiles_target-pheno.phen

#Set Phenotype column in .fam file to -9
mv $output_dir$prefix_bfiles_target-target-qced.fam $output_dir$prefix_bfiles_target-target-qced.fam.old
awk '{ $6 = -9; print }' $output_dir$prefix_bfiles_target-target-qced.fam.old > $output_dir$prefix_bfiles_target-target-qced.fam


## 2. Sex covariate
#Extract Sex covariate from .fam file and combine with PCA file
awk '{print $1, $2, $5}' $output_dir$prefix_bfiles_target-target-qced.fam > $output_dir$prefix_bfiles_target-qced-sex-temp.covar
awk 'BEGIN {print "FID IID Sex"} 1' $output_dir$prefix_bfiles_target-qced-sex-temp.covar > $output_dir$prefix_bfiles_target-qced-sex.covar
rm $output_dir$prefix_bfiles_target-qced-sex-temp.covar

## Concatenate Covariates into one file -> PCA and Sex covariate
# Add Header to files
awk 'BEGIN {print "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10"} {print $0}' \
     $output_dir$prefix_bfiles_target-grm-pca.eigenvec > temp_pca.covar
Rscript $SCRIPTS_DIR/merge_covariates.R temp_pca.covar $output_dir$prefix_bfiles_target-qced-sex.covar $output_dir$prefix_bfiles_target-qced.covariate
rm temp_pca.covar

