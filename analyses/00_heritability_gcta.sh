#!/bin/bash
# This script filters and prunes genotype data, then calculates principal components (PCA) and a genetic relationship matrix (GRM) using PLINK and GCTA.
# It then estimates heritability for the full dataset, as well as separately for SNPs within lead regions and SNPs outside those lead regions.

## Setting up variables 
geno_dir='/project/geno/'
output_dir='/project/output/'
pheno_dir='/project/pheno/'
scripts_dir='/project/scripts/'
leads_region_dir='/project/h2-leads-region/'
nonleads_region_dir='/project/h2-nonleads-region/'
prefix_bfiles='filename_prefix'
LDregions_file='high-LD-regions.txt' 
prefix_outfile_grm='outfile_name'
pheno_file='pheno_filename'
leads_regionSNPs_file='keep_SNPs_lead_region.txt'
nonleads_regionSNPs_file='keep_SNPs_outof_lead_region.txt'
threadnum=50
memory_num=300000

## Get Phenotype IDS
awk '{print $1, $2}' $pheno_dir$pheno_file > $output_dir$pheno_file-ids

## 1.Drop samples without phenotypes from  bfiles
plink2 --bfile $geno_dir$prefix_bfiles \
       --keep $output_dir$pheno_file-ids \
       --make-bed \
       --out $output_dir$prefix_bfiles

# Check how many cases and controls
echo "Number of cases and controls"
awk '{print $6}' $geno_dir$prefix_bfiles.fam | sort | uniq -c

## Sex covariate
cut -f1,2,5 $geno_dir$prefix_bfiles.fam > $output_dir${prefix_outfile_grm}-sex.covar


## 2. Prune Geno DATA for computing GRM and PCA
## a. Excluding high-LD and HLA regions
plink --bfile $output_dir$prefix_bfiles \
      --chr 1-22 \
      --maf 0.01 \
      --geno 0.1 \
      --hwe 5e-6 \
      --exclude $LDregions_file \
      --indep-pairwise 50 5 0.2 \
      --memory 200000 \
      --threads ${threadnum} \
      --make-bed \
      --out $output_dir$prefix_bfiles-pruning \

## b. LD prune
plink --bfile $output_dir$prefix_bfiles-pruning \
      --extract $output_dir$prefix_bfiles-pruning.prune.in \
      --make-bed \
      --out $output_dir$prefix_bfiles-LDpruned \
      --memory ${memory_num} \
      --threads ${threadnum}

## 3. Generate PCA with prunned data
plink --bfile $output_dir$prefix_bfiles-LDpruned \
      --pca 10 \
      --memory ${memory_num} \
      --threads ${threadnum} \
      --out $output_dir$prefix_outfile_grm-pca
     

## 4. Compute GRM with prunned data
gcta64 --bfile $output_dir$prefix_bfiles-LDpruned \
       --make-grm --thread-num 10 \
       --out $output_dir$prefix_outfile_grm-grm


#### (I) FULL DATA Heritability

## Heritability Early Age at Onset (< 50 years)
gcta64 --reml \
       --grm $output_dir$prefix_outfile_grm-grm \
       --pheno $pheno_dir$pheno_file \
       --prevalence 0.002 \
       --qcovar $output_dir$prefix_outfile_grm-pca.eigenvec \
       --covar $output_dir$prefix_outfile_grm-sex.covar\
       --out $output_dir$prefix_outfile_grm-heritability-early-onset \
       --thread-num 200



### (II) LEAD SNP REGIONS Heritability OR PD GWAS SNPs in PD loci/regions
# -> Extract and use only SNPs within +/- 1Mb region of independent lead SNPs
# -> Re-compute GRM and estimate heritability
# -> Use PCA computed from all genotypes

## a. Extract Geno Data in Lead Regions
plink2 --bfile $output_dir$prefix_bfiles \
       --extract $geno_dir$leads_regionSNPs_file \
       --make-bed \
       --out $leads_region_dir$prefix_bfiles

## b. Compute GRM with extracted Sequence data
gcta64 --bfile $leads_region_dir$prefix_bfiles \
       --make-grm --thread-num 10 \
       --out $leads_region_dir$prefix_outfile_grm-grm

## c. Compute Heritability
gcta64 --reml \
       --grm $leads_region_dir$prefix_outfile_grm-grm \
       --pheno $pheno_dir$pheno_file \
       --prevalence 0.002 \
       --qcovar $output_dir$prefix_outfile_grm-pca.eigenvec \
       --covar $output_dir$prefix_outfile_grm-sex.covar\
       --out $leads_region_dir$prefix_outfile_grm-heritability-early-onset-PD-GWAS-loci \
       --thread-num 200



#### (III) NON LEAD SNP REGIONS Heritability OR PD GWAS Excluded
# -> Extract and use all SNPs not within +/- 1Mb region of independent lead SNPs 
# -> Re-Compute GRM and estimate heritability

## a. Extract Geno Data in Lead Regions
plink2 --bfile $output_dir$prefix_bfiles \
       --extract $geno_dir$nonleads_regionSNPs_file \
       --make-bed \
       --out $nonleads_region_dir$prefix_bfiles

## b. Compute GRM with extracted Sequence data
gcta64 --bfile $nonleads_region_dir$prefix_bfiles \
       --make-grm --thread-num 10 \
       --out $nonleads_region_dir$prefix_outfile_grm-grm

## c. Compute Heritability
gcta64 --reml \
       --grm $nonleads_region_dir$prefix_outfile_grm-grm \
       --pheno $pheno_dir$pheno_file \
       --prevalence 0.002 \
       --qcovar $output_dir$prefix_outfile_grm-pca.eigenvec \
       --covar $output_dir$prefix_outfile_grm-sex.covar\
       --out $nonleads_region_dir$prefix_outfile_grm-heritability-early-onset-PD-GWAS-loci-excluded \
       --thread-num 200



