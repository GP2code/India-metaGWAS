
#!/usr/bin/env Rscript

# Script: Differential Gene Expression Analysis
#
# Description:
#This script calculates p-adjusted values for pairwise comparisons of gene expression
# medians between superpopulations using a permutation test method. It computes the median
# expression for each gene within each superpopulation and evaluates the significance of 
# differences in these medians across superpopulations.

# Required Input Files:
# 1. Expression Data (BED format): Contains DESeq2 normalized gene expression values from the previous study "Taylor, D.J., Chhetri, S.B., Tassia, M.G. et al. Sources of gene expression variation in a globally diverse human cohort. Nature 632, 122–130 (2024). https://doi.org/10.1038/s41586-024-07708-2".
# 2. Annotation File (TXT format): Includes metadata for sample populations (Taylor, D.J., Chhetri, S.B., Tassia, M.G. et al. Sources of gene expression variation in a globally diverse human cohort. Nature 632, 122–130 (2024). https://doi.org/10.1038/s41586-024-07708-2).
# 3. Gene Info File (CSV format): Maps Ensembl IDs to gene symbols and SNP information.
#
# Outputs:
# - A CSV file containing p-values, gene medians, and metadata for each comparison.
#
# Usage:
# Update the file paths under the "Input files" section before running the script.


expression_file <- "path/to/expression_data.bed"
annotation_file <- "path/to/sample_metadata.txt"
gene_info_file <- "path/to/gene_info.csv"
output_file <- "path/to/permute_test_result.csv"

# Error handling for input files
if (!file.exists(expression_file)) stop("Expression data file not found. Please check the file path.")
if (!file.exists(annotation_file)) stop("Annotation file not found. Please check the file path.")
if (!file.exists(gene_info_file)) stop("Gene info file not found. Please check the file path.")



library(ggplot2)
library(dplyr)
library(tidyr)

# Set up permutation test functions
perm_test_median <- function(group1, group2, n_perm = 100000) {
  
  combined <- c(group1, group2)  
  n1 <- length(group1)
  n2 <- length(group2)
  
  obs_diff <- median(group1) - median(group2)
  
  perm_diffs <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    permuted <- sample(combined)  # Randomly shuffle the combined data
    perm_group1 <- permuted[1:n1]  # First n1 elements go to group 1
    perm_group2 <- permuted[(n1+1):(n1+n2)]  # Rest go to group 2
    perm_diffs[i] <- median(perm_group1) - median(perm_group2)
  }
  # Calculate p-value
  p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
  
  if (p_value == 0) {
    p_value <- 1 / n_perm
  }
  return(p_value)
} 


# Load data
expression_data <- read.table(expression_file, header = TRUE, sep = "\t")
annotation <- read.table(annotation_file, header = TRUE, sep = "\t")
gene_info <- read.csv(gene_info_file)


# Prepare count matrix
counts_matrix <- expression_data[,5:ncol(expression_data)]
rownames(counts_matrix) <- expression_data[,4]

# Transpose the expression matrix
expression_matrix_t <- as.data.frame(t(counts_matrix))
expression_matrix_t$SampleID <- rownames(expression_matrix_t)

# Prepare Annotation File
annotation <- data.frame(annotation$sample_kgpID, annotation$continentalGroup, annotation$population)
colnames(annotation) <- c("SampleID", "SuperPopulation", "Population")
annotation$SuperPopulation[which(annotation$Population == "FIN")] <- "FIN"

# Merge expression data with annotation
merged_data <- merge(expression_matrix_t, annotation, by = "SampleID")

# Convert the data to long format
long_data <- merged_data %>%
  pivot_longer(cols = -c(SampleID, Population, SuperPopulation), 
               names_to = "Gene", values_to = "Expression")
long_data$Gene <- gsub("\\..*", "", long_data$Gene)


# Consolidate all results in a single data frame
all_results <- data.frame(SNP = character(), Gene = character(),
                          SuperPopulation1 = character(), SuperPopulation2 = character(),
                          P_Value = numeric(), stringsAsFactors = FALSE)

# Permutation test for each gene
for (i in 1:nrow(gene_info)) {
  ensembl_id <- gene_info$EnsemblID[i]
  gene_symbol <- gene_info$GeneSymbol[i]
  snp_id <- gene_info$SNP[i]  
  
  gene_data <- long_data %>% filter(Gene == ensembl_id)
  
  # Check if there are any samples for this gene
  if (nrow(gene_data) > 0) {
    groups <- unique(gene_data$SuperPopulation)
    
    # Perform pairwise permutation tests between all pairs of super populations
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        group1_values <- gene_data$Expression[gene_data$SuperPopulation == groups[i]]
        group2_values <- gene_data$Expression[gene_data$SuperPopulation == groups[j]]
        
        # Perform permutation test
        p_value <- perm_test_median(group1_values, group2_values)
        
        
        all_results <- rbind(all_results, data.frame(SNP = snp_id,
                                                     Gene = gene_symbol, 
                                                     SuperPopulation1 = groups[i], 
                                                     SuperPopulation2 = groups[j], 
                                                     P_Value = p_value))
      }
    }
  }
}

# Adjust p-values using the Benjamini-Hochberg procedure
all_results <- all_results %>%
  group_by(Gene) %>%
  mutate(P_adjusted = p.adjust(`P_Value`, method = "fdr"))

# Add median expression per gene for each superpopulation

gene_medians <- long_data %>%
  group_by(Gene, SuperPopulation) %>%
  summarize(Median_Expression = median(Expression, na.rm = TRUE), .groups = "drop")

gene_medians <- gene_medians %>%
  left_join(gene_info, by = c("Gene" = "EnsemblID")) %>%
  select(GeneSymbol, SuperPopulation, Median_Expression) %>%
  rename(Gene = GeneSymbol)

final_results <- all_results %>%
  left_join(gene_medians, by = c("Gene", "SuperPopulation1" = "SuperPopulation")) %>%
  rename(Median_Expression_SP1 = Median_Expression) %>%
  left_join(gene_medians, by = c("Gene", "SuperPopulation2" = "SuperPopulation")) %>%
  rename(Median_Expression_SP2 = Median_Expression)

write.csv(final_results, output_file, row.names = FALSE)