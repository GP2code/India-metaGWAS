#!/usr/bin/env Rscript

# This script generates plots based on p-adjust values and gene expression.

# Required Input Files:
# 1. Expression Data (BED format): Contains DESeq2 normalized gene expression values from the previous study "Taylor, D.J., Chhetri, S.B., Tassia, M.G. et al. Sources of gene expression variation in a globally diverse human cohort. Nature 632, 122–130 (2024). https://doi.org/10.1038/s41586-024-07708-2".
# 2. Annotation File (TXT format): Includes metadata for sample populations (Taylor, D.J., Chhetri, S.B., Tassia, M.G. et al. Sources of gene expression variation in a globally diverse human cohort. Nature 632, 122–130 (2024). https://doi.org/10.1038/s41586-024-07708-2).
# 3. Gene Info File (CSV format): Maps Ensembl IDs to gene symbols and SNP information.
# 4. significance_data file: the output file (permute_test_result.csv) from the Differential_Gene_Expression_Analysis code

# Outputs:
# - Saves plots in a PDF file.
#
# Usage:
# Update the file paths under the "Input files" section before running the script.


expression_file <- "path/to/expression_data.bed"
annotation_file <- "path/to/sample_metadata.txt"
gene_info_file <- "path/to/gene_info.csv"
significance_data <- "path/to/permute_test_result"
output_file <- "path/to/violin_boxPlot.pdf"

# Error handling for input files
if (!file.exists(expression_file)) stop("Expression data file not found. Please check the file path.")
if (!file.exists(annotation_file)) stop("Annotation file not found. Please check the file path.")
if (!file.exists(gene_info_file)) stop("Gene info file not found. Please check the file path.")
if (!file.exists(significance_data)) stop("significance_data file not found. Please check the file path.")




library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(gridExtra) 
library(cowplot) 
library(patchwork)

# Prepare count matrix
counts_matrix <- expression_file[,5:ncol(expression_file)]
rownames(counts_matrix) <- expression_file[,4]

# Transpose expression matrix to merge with annotation
expression_matrix_t <- as.data.frame(t(counts_matrix))
expression_matrix_t$SampleID <- rownames(expression_matrix_t)

# prepare annotation data
annotation <- data.frame(annotation$sample_kgpID, annotation$continentalGroup, annotation$population)
colnames(annotation) <- c("SampleID", "SuperPopulation", "Population")
annotation$SuperPopulation[annotation$Population == "FIN"] <- "FIN"

# Merge expression data with annotation
merged_data <- merge(expression_matrix_t, annotation, by = "SampleID")

# Convert data to long format
long_data <- merged_data %>%
  pivot_longer(cols = -c(SampleID, Population, SuperPopulation), 
               names_to = "Gene", values_to = "Expression")
long_data$Gene <- gsub("\\..*", "", long_data$Gene)

# Define population labels
pop_labels <- long_data %>%
  select(Population, SuperPopulation) %>%
  distinct() %>%
  mutate(Label = paste0(Population, "\n(", SuperPopulation, ")"))

pop_labels <- pop_labels %>%
  arrange(SuperPopulation, Population)

# Initialize a list to store all plots
plot_list <- list()

# Loop through the genes
for (i in 1:nrow(gene_info)) {  # Change to your desired number of genes or use nrow(gene_info)
  ensembl_id <- gene_info$EnsemblID[i]
  gene_symbol <- gene_info$GeneSymbol[i]
  snp_id <- gene_info$SNP[i]  # Get SNP ID for the title
  
  # Filter data for the current gene
  gene_data <- long_data %>% filter(Gene == ensembl_id)
  
  # Only proceed if gene_data is not empty
  if (nrow(gene_data) > 0) {
    gene_data$Population <- factor(gene_data$Population, levels = pop_labels$Population)
    
    # Create the plot 
    p <- ggplot(gene_data, aes(x = SuperPopulation, y = Expression, fill = SuperPopulation)) +
      geom_violin(alpha = 0.6, trim = FALSE, show.legend = FALSE) +  # Hide legend in this layer
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, show.legend = FALSE) +  # Hide legend in this layer
      ggtitle(label = snp_id,subtitle = gene_symbol) +  # Add SNP ID as title and Gene Symbol as subtitle
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",  # Remove legend from each individual plot
        plot.title = element_text(hjust = 0.5, size = 14,face = "bold"),  # Title in the top center (bold SNP ID)
        plot.subtitle = element_text(hjust = 0, size = 12)  # Subtitle in the top left (Gene Symbol)
      ) +
      labs(
        x = "Super Population",
        y = "Expression Count"
      )
    
    # Add significance layers if available
    gene_significance <- significance_data %>% filter(Gene == gene_symbol)
    if (nrow(gene_significance) > 0) {
      base_y <- max(gene_data$Expression) + 0.5
      step <- 0.3  # Distance between each layer
      
      for (j in 1:nrow(gene_significance)) {
        pop1 <- gene_significance$SuperPopulation1[j]
        pop2 <- gene_significance$SuperPopulation2[j]
        P_Value <- gene_significance$P_adjusted[j]
        
        annotation_text <- ifelse(P_Value < 0.001, "***",
                                  ifelse(P_Value < 0.01, "**",
                                         ifelse(P_Value < 0.05, "*", "")))
        
        if (annotation_text != "") {
          y_position <- base_y + (j - 1) * step
          p <- p + geom_signif(
            comparisons = list(c(pop1, pop2)),
            annotations = annotation_text,
            y_position = y_position,
            size = 0.3,
            textsize = 5,
            color = "black"
          )
        }
      }
    }
    
    # Append the valid plot to plot_list
    plot_list <- append(plot_list, list(p))
  }
}


number_column=3 # this is based on number of column in pdf file; adjusting it based it based on the number of plots.

# using patchwork to combine the plots into a 3x3 grid with a shared legend
combined_plot <- wrap_plots(plotlist = plot_list, ncol = number_column) + 
  theme(legend.position = "none")  # Remove legend from individual plots


pdf(output_file, width = 15, height = 15)  # Adjusting width and height based on the number of plots
print(combined_plot)
dev.off()
