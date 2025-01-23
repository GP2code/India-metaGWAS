# Deciphering the Genetic Architecture of Parkinson’s Disease in India

`GP2 ❤️ Open Science 😍`

Pending DOI
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


**Last Updated:** January 2025

## Summary

This is the online repository for the manuscript titled **"Deciphering the Genetic Architecture of Parkinson’s Disease in India"**. This study represents the first largest genome-wide assessment of Parkinson’s disease in the Indian population.

## Data Statement


## Workflow

1. Quality cotrol (QC) and genotype processing with IDEAL-GENOM-QC pipeline: https://github.com/cge-tubingens/IDEAL-GENOM-QC
2. Imputation to whole genome sequence with the Michigan imputation server
3. GWAS with IDEAL-GENOM pipeline: https://github.com/cge-tubingens/IDEAL-GENOM
4. Heritability estimation with GCTA
5. Meta-Analysis with GWAMA
6. MAPT haplotype analysis with haplo.stats
7. Gene and tissue enrichment analyses with FUMA/MAGMA
8. Polygenic risk prediction with PRSice-2
9. Gene expression analyses
    

## Repository Orientation
* The `analyses/` directory includes all post GWAS analyses, for QC+GWAS, see the IDEAL-GENOM-QC pipeline: https://github.com/cge-tubingens/IDEAL-GENOM-QC and https://github.com/cge-tubingens/IDEAL-GENOM

```
analyses/
├── 00_heritability_gcta.sh
├── 01_harmonize_meta_sumstats_popAB.py
├── 02_meta_data_prep_GWAMA.R
├── 03_meta_analysis_GWAMA.sh
├── 04_meta_COJO_gcta.R
├── 05_prs_BaseDataQC_IndianCohort.sh
├── 06_prs_TargetDataQC.sh
├── 07_prs_analysis_indianBase.sh
├── 08_prs_ROC_AUC_plots.R
├── 09_mapt_haplotype_analysis.ipynb
├── 10_gene_expression_analysis.R
└── 11_gene_expression_plots.R
```

## Analysis Scripts & Notebooks
* Languages: Python, R, bash

| **Notebooks  / Scripts**            | **Description**                                                                            |
|:-----------------------------------:|:------------------------------------------------------------------------------------------:|
| `00_heritability_gcta.sh`             | Estimating heritability with SNP data                                                      |
| `01_harmonize_meta_sumstats_popAB.py` | Harmonizing GWAS summary statistics for meta-analyses                                      |
| `02_meta_data_prep_GWAMA.R`           | Standardizes two GWAS summary stats for meta-analysis in GWAMA                            |
| `03_meta_analysis_GWAMA.sh`           | Running meta-analysis with Indian GWAS and Multi-ethic  GWAS                               |
| `04_meta_COJO_gcta.R`                 | Conditional association analysis on meta output                                             |
| `05_prs_BaseDataQC_IndianCohort.sh`   | Quality control of base data for PRS analyses                                              |
| `06_prs_TargetDataQC.sh`              | Quality control of target data for PRS analyses                                            |
| `07_prs_analysis_indianBase.sh`       | Conducting PRS analysis with  Luxgiant variants and summary statistics Andrews et al.,2023 |
| `08_prs_ROC_AUC_plots.R`              | ROC (Receiver Operating Characteristic) and AUC (Area Under the Curve) for PRS             |
| `09_mapt_haplotype_analysis.ipynb`    | *MAPT* loci haplotype analyses                                                               |
| `10_gene_expression_analysis.R`       | Differential gene expression analyses                                                      |
| `11_gene_expression_plots.R`          | Visualization gene expression analyses                                                     |



## Software

| **Software**                        | **Version(s)**  | **Resource URL**                                                   | **RRID**        | **Notes**                                                                                                                                                     |
|:-----------------------------------:|:---------------:|:------------------------------------------------------------------:|:---------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| Python Programming Language         | 3.10 and 3.12   | https://www.python.org/                                            | RRID:SCR_008394 | ideal-genom, ideal-genom-qc; libraries that wraps quality control, visualization, and GWAS                                                                    |
| PLINK                               | 1.9 and 2.0     | http://www.nitrc.org/projects/plink                                | RRID:SCR_001757 | GWAS with GLM                                                                                                                                                 |
| samtools (bcftools)                 | 1.2             | https://www.htslib.org/                                            | RRID:SCR_002105 | VCF manipulation                                                                                                                                              |
| GCTA                                | 1.94.1          | https://yanglab.westlake.edu.cn/software/gcta/#Overview            | not available   | GWAS with mixed model, estimating heritability                                                                                                                |
| PRSice-2                            | 2.3.5           | https://choishingwan.github.io/PRSice/                             | not available   | used to perform PRS analyses                                                                                                                                  |
| GWAMA                               | 2.2.2           | https://genomics.ut.ee/en/tools                                    | RRID:SCR_006624 | used to perform GWAS meta analyses                                                                                                                            |
| Michigan Imputation Server          | 2.0.0           | https://imputationserver.readthedocs.io/en/latest/getting-started/ | RRID:SCR_023554 | used to impute genotype data to whole genome sequence                                                                                                         |
| R Project for Statistical Computing | 4.3.3 and 4.4.1 | http://www.r-project.org/                                          | RRID:SCR_001905 | pROC, haplo.stats, tidyr, ggplot2, dplyr, ggsignif, gridExtra, cowplot, patchwork; visualization,  analysis of indirectly measured haplotypes, data wrangling |
