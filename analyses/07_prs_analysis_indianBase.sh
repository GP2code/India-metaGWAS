#!/bin/bash
## PRS Luxgiant (Target) - Indian GWAS (Base)
# Conducting PRS analysis with LuxGiant variants and summary statistics Andrews et al.,2023

## Directories
PRSCIE_EXEC='/home/ubuntu/bin/PRSice/'
inputdir_basedata='/project/indian-cohort/'
inputdir_targetdata='/project/target/'
working_dir='/project/..'
output_one_dir='/project/run-one/'
output_two_dir='/project/run-two/'
output_three_dir='/project/run-three/'
prefix_bfiles_target='filename_prefix'
basefile='filename.txt'
name_prsfile='filenae.txt'


#=========================================================
### Run One 
## Best  performing model from thresholding and clumping
#=========================================================
## With lowest quantile as reference
Rscript $PRSCIE_EXEC/PRSice.R \
    --prsice $PRSCIE_EXEC/PRSice_linux \
    --base $inputdir_basedata$basefile \
    --target $inputdir_targetdata$prefix_bfiles_target-target-qced \
    --clump-kb 500kb \
    --clump-r2 0.500000 \
    --pheno $inputdir_targetdata$prefix_bfiles_target-pheno.phen \
    --cov $inputdir_targetdata$prefix_bfiles_target-qced.covariate \
    --pheno-col PD \
    --binary-target T \
    --prevalence 0.005 \
    --stat OR \
    --or \
    --quantile 4 \
    --quant-ref 1 \
    --print-snp  \
    --seed 12345 \
    --thread 50 \
    --out $output_one_dir${name_prsfile}


## Run Two
# Run with permutations
Rscript $PRSCIE_EXEC/PRSice.R \
    --prsice $PRSCIE_EXEC/PRSice_linux \
    --base $inputdir_basedata$basefile \
    --target $inputdir_targetdata$prefix_bfiles_target-target-qced \
    --clump-kb 500kb \
    --clump-r2 0.500000 \
    --pheno $inputdir_targetdata$prefix_bfiles_target-pheno.phen \
    --cov $inputdir_targetdata$prefix_bfiles_target-qced.covariate \
    --binary-target T \
    --prevalence 0.005 \
    --stat OR \
    --or \
    --quantile 5 \
    --quant-ref 1 \
    --bar-levels 5e-06,0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    --logit-perm T \
    --perm 10000 \
    --thread 50 \
    --out $output_two_dir${name_prsfile}


#=========================================================
### Run Three
## 2nd bset performing model from thresholding and clumping
#=========================================================
## With lowest quantile as reference
Rscript $PRSCIE_EXEC/PRSice.R \
    --prsice $PRSCIE_EXEC/PRSice_linux \
    --base $inputdir_basedata$basefile \
    --target $inputdir_targetdata$prefix_bfiles_target-target-qced \
    --clump-kb 500kb \
    --clump-r2 0.500000 \
    --pheno $inputdir_targetdata$prefix_bfiles_target-pheno.phen \
    --cov $inputdir_targetdata$prefix_bfiles_target-qced.covariate \
    --binary-target T \
    --prevalence 0.005 \
    --stat OR \
    --or \
    --quantile 4 \
    --quant-ref 1 \
    --bar-levels 5e-06 \
    --fastscore \
    --all-score \
    --print-snp  \
    --seed 12345 \
    --thread 50 \
    --out $output_three_dir${name_prsfile}



