#!/usr/bin/env Rscript
# ROC (Receiver Operating Characteristic) and AUC (Area Under the Curve) for PRS

compute_and_plot_roc_test_sig <- function(pheno_data, prs_data, covar_data, plot_filename = "roc_curve.jpeg", results_filename = "roc_auc_results.tsv",
                                 ancestry=NULL, get_roc_prs=FALSE) {
 #Function to compute ROC, AUC, plot figures of Full and Null models annotated with CI, and compute significance
  library(pROC)

  # Step 1: Merge the phenotype, PRS, and covariates data using common columns (e.g., FID, IID)
  merged_data <- merge(pheno_data, prs_data, by = c("FID", "IID"))
  merged_data <- merge(merged_data, covar_data, by = c("FID", "IID"))

  if (get_roc_prs) {
    # Full model: PRS + sex + 10 PCs
    full_model <- glm(PD ~ PRS + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      data = merged_data, family = binomial)
    
    # Null model: covariates only (without PRS)
    null_model <- glm(PD ~ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      data = merged_data, family = binomial)
    
    # PRS only model
    model_prs <- glm(PD ~ PRS, data = merged_data, family = binomial)

    # Compute ROC curves
    roc_full <- roc(merged_data$PD, fitted(full_model))
    roc_null <- roc(merged_data$PD, fitted(null_model))
    roc_prs <- roc(merged_data$PD, fitted(model_prs))

    # Compute AUCs and Confidence Intervals
    auc_full <- auc(roc_full)
    auc_null <- auc(roc_null)
    auc_prs <- auc(roc_prs)

    ci_full <- ci.auc(roc_full)
    ci_null <- ci.auc(roc_null)
    ci_prs <- ci.auc(roc_prs)

    # Compute p-values for the significance of ROC differences
    p_value_full_vs_null <- roc.test(roc_full, roc_null)$p.value  
    p_value_full_vs_prs <- roc.test(roc_full, roc_prs)$p.value

    # Plot ROC curves
    jpeg(plot_filename, width = 800, height = 800)
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, cex.sub = 1.5)
    plot(roc_full, col = "red", lwd = 3, main = "")
    plot(roc_prs, col = "turquoise", add = TRUE, lwd = 3)
    plot(roc_null, col = "orange", add = TRUE, lwd = 3)
    legend("topleft", 
           legend = c(sprintf("AUC (PRS + COV) = %.2f (95%% CI: %.2f - %.2f)", auc_full, ci_full[1], ci_full[3]),
                      sprintf("AUC (PRS) = %.2f (95%% CI: %.2f - %.2f)", auc_prs, ci_prs[1], ci_prs[3]),
                      sprintf("AUC (COV) = %.2f (95%% CI: %.2f - %.2f)", auc_null, ci_null[1], ci_null[3])),
           col = c("red", "turquoise", "orange"), lwd = 3, cex = 1.2)
    dev.off()

    # Collect ROC, AUC, and p-values
    results <- data.frame(
      Model = c("PRS + COV", "PRS", "COV"),
      AUC = c(auc_full, auc_prs, auc_null),
      CI_Lower = c(ci_full[1], ci_prs[1], ci_null[1]),
      CI_Upper = c(ci_full[3], ci_prs[3], ci_null[3]),
      P_Value_Comparison = c(p_value_full_vs_null, p_value_full_vs_prs, NA) # p-value for each model comparison
    )
    results$Ancestry <- ancestry

  } else {
    # Full model: PRS + sex + 10 PCs
    full_model <- glm(PD ~ PRS + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      data = merged_data, family = binomial)
    # Null model: covariates only (without PRS)
    null_model <- glm(PD ~ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      data = merged_data, family = binomial)

    # Compute ROC curves
    roc_full <- roc(merged_data$PD, fitted(full_model))
    roc_null <- roc(merged_data$PD, fitted(null_model))

    # Compute AUCs and Confidence Intervals
    auc_full <- auc(roc_full)
    auc_null <- auc(roc_null)

    ci_full <- ci.auc(roc_full)
    ci_null <- ci.auc(roc_null)

    # Compute p-value for significance of ROC difference
    p_value_full_vs_null <- roc.test(roc_full, roc_null)$p.value

    # Plot ROC curves
    jpeg(plot_filename, width = 800, height = 800)
    par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, cex.sub = 1.5)
    plot(roc_full, col = "red", lwd = 3, main = "")
    plot(roc_null, col = "turquoise", add = TRUE, lwd = 3)
    legend("topleft", 
           legend = c(sprintf("AUC (PRS + COV) = %.2f (95%% CI: %.2f - %.2f)", auc_full, ci_full[1], ci_full[3]),
                      sprintf("AUC (COV) = %.2f (95%% CI: %.2f - %.2f)", auc_null, ci_null[1], ci_null[3])),
           col = c("red", "turquoise"), lwd = 3, cex = 1.2)
    dev.off()

    # Collect ROC, AUC, and p-value
    results <- data.frame(
      Model = c("Full Model", "Null Model"),
      AUC = c(auc_full, auc_null),
      CI_Lower = c(ci_full[1], ci_null[1]),
      CI_Upper = c(ci_full[3], ci_null[3]),
      P_Value_Comparison = c(p_value_full_vs_null, NA) # p-value only for the comparison
    )
    results$Ancestry <- ancestry
  }

  # Save results to a specified text file
  write.table(results, file = results_filename, sep = "\t", row.names = FALSE, quote = FALSE)

  # Print AUC values and p-values for reference
  print(paste("AUC (Full Model):", auc_full))
  print(paste("AUC (Null Model):", auc_null))
  print(paste("P-value for AUC comparison (Full vs Null):", p_value_full_vs_null))
  if (get_roc_prs) print(paste("P-value for AUC comparison (Full vs PRS):", p_value_full_vs_prs))
}
#==================
# Example usage:
#=====================
#compute_and_plot_roc_test_sig(pheno_data, prs_data=prs_dat_eu, covar_data, plot_filename = "test_roc_plot_european_significance.jpeg",
#                              results_filename = "test2_table_roc_results_european_significance.tsv", ancestry='EUR', get_roc_prs=TRUE)

#compute_and_plot_roc_test_sig(pheno_data, prs_data=prs_dat2_ind, covar_data, plot_filename = "test_roc_plot_indian_significance.jpeg",
#                              results_filename = "test2_table_roc_results_indian_significance.tsv", ancestry='IND', get_roc_prs=TRUE)


# Function to plot Odds Ratios for Quartiles with Confidence Intervals
plot_quartiles <- function(data, output_filename, num_snps1='', num_snps2='',filename_dataframe_out = NULL) {
  # Exclude the reference group
  data_filtered <- data[data$Quantile != 1, ]  # Exclude reference 

  # Plot Odds Ratios with Confidence Intervals
  plot <- ggplot(data_filtered, aes(x = factor(Quantile), y = OR, color = Ancestry)) +
    geom_point(size = 3) +  # Plot the OR as points
    geom_errorbar(aes(ymin = CI.L, ymax = CI.U), width = 0.2) +  # Add error bars for CIs
    labs(title = "",
         x = "Quartile of PRS",   # Custom x-axis label
        y = expression(paste("Odds ratios (95% CI), log"[10], " scaled")), color = "Base set") +  
    scale_x_discrete(labels = c("2nd", "3rd", "4th")) +  # Custom labels for x-axis
    scale_y_log10() +  # Use log scale for Odds Ratio
    theme_minimal(base_size = 15) +  # Minimal theme
    theme(panel.grid = element_blank(),  # Remove all grid lines
          panel.background = element_rect(fill = "white"),  # Set plain white background
          plot.background = element_rect(fill = "white"),  # Set plot background to white
          panel.border = element_blank(),  # Remove outer border frame
          axis.line = element_line(),    # Keep the axis lines
          axis.ticks = element_line(), # Keep the axis ticks
          legend.position = c(0.1, 1),  # Move legend to top-left
          legend.justification = c("left", "top")) +  # Align legend to top-left corner
    scale_color_manual(values = c("EUR" = "turquoise", "IND" = "coral"))  # Custom colors for each ancestry group

  # Save the plot 
  ggsave(output_filename, plot, width = 7, height = 7, dpi = 300)

  # Optionally save the input data as a tab-delimited file
  if (!is.null(filename_dataframe_out)) {
    write.table(data, file = filename_dataframe_out, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

#==========================================
# Example of how to use the function
#=========================================
#plot_quartiles(data, "test_odds_ratios_quantiles_prs_european.jpeg", 'test_table_odds_ratios_data_european.tsv')
#plot_odds_ratios(data, "odds_ratios_quantiles_prs_european.jpeg", "odds_ratios_data.tsv")


plot_prs_distribution <- function(data, plot_filename = "prs_distribution.jpeg", save.plot=NULL) {
  ggplot(data, aes(x = PRS, fill = factor(PD))) +
    geom_density(alpha = 0.5) +
    labs(title = "",
         x = "Polygenic Risk Score (PRS)",
         y = "Density",
         fill = "Phenotype") +
    scale_fill_manual(values = c("blue2", "red"), labels = c("Control", "Case")) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",          # Position the legend at the top
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.text = element_text(size = 16),  # Increase axis tick labels size
    axis.title = element_text(size = 16)  # Increase axis labels size
  )
  
  # Save the plot
  if(!is.null(save.plot)){
      ggsave(plot_filename, width = 10, height = 7)
  }
  
}

##---------------------------------------------------
plot_prs_distribution <- function(data, plot_filename="prs_distribution_plot.jpeg", save.plot=NULL) {
  # Ensure the output file has a .jpeg extension
  if (tools::file_ext(plot_filename) != "jpeg") {
    stop("Output file must have a .jpeg extension.")
  }
  
  # Create the box plot
  box_plot <- ggplot(data, aes(x = factor(PD), y = PRS, fill = factor(PD))) +
    geom_boxplot(alpha = 0.7) +
    labs(x = "", y = "PRS") +
    scale_fill_manual(values = c("royalblue3", "orange3"), labels = c("Control", "Case")) +
    scale_x_discrete(labels = c("Control", "Case")) +  # Change x-axis tick labels
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "none",           # Remove legend
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.text = element_text(size = 16),  # Increase axis tick labels size
      axis.title = element_text(size = 16)  # Increase axis labels size
    )
  
  # Save the plot
  if(!is.null(save.plot)){
      ggsave(plot_filename, width = 10, height = 7)
  }
  
}

# Example usage:
# Assuming prs_pheno_eu is the input data frame
#plot_prs_distribution(prs_pheno_eu, "prs_distribution_plot.jpeg")