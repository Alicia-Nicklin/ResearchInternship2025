# Method comparison - z-score vs ssGSEA 
# Task 3 from my internship - comparing the two approaches Daniel showed me
# Using base R plots because ggsurvplot kept breaking

library(survival)
library(dplyr)
library(ggplot2)

# Load my data as usual
setwd("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA")
DF <- read.csv("merged_COAD_primary_4genes_vFINAL.csv")

# Make folders to organize everything - learned this is important!
main_folder <- "KM_task3_method_comparison"
zscore_folder <- file.path(main_folder, "zscore_plots")
ssgsea_folder <- file.path(main_folder, "ssgsea_plots")
comparison_folder <- file.path(main_folder, "comparison_results")

for(folder in c(main_folder, zscore_folder, ssgsea_folder, comparison_folder)) {
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
    print(paste("Created folder:", folder))
  }
}

# Settings - using median split like we discussed
THRESHOLD <- 0.50
gene_signatures <- list(
  "CES1_CPT1A" = c("CES1", "CPT1A"),
  "CES1_MGLL_ACSL3" = c("CES1", "MGLL", "ACSL3"), 
  "CES1_MGLL_CPT1A" = c("CES1", "MGLL", "CPT1A"),
  "CES1_CPT1A_ACSL3" = c("CES1", "CPT1A", "ACSL3"),
  "CES1_MGLL_CPT1A_ACSL3" = c("CES1", "MGLL", "CPT1A", "ACSL3")
)

print("Method comparison analysis")
print(paste("Working with", nrow(DF), "patients"))

# Method 1: Z-score approach (using the columns I already calculated)
get_zscore_group <- function(df, sig_name, threshold = 0.50) {
  zscore_col <- paste0("zsign_", sig_name)
  if(!zscore_col %in% colnames(df)) {
    print(paste("Can't find z-score column:", zscore_col))
    return(NULL)
  }
  
  zscore_values <- df[[zscore_col]]
  threshold_value <- quantile(zscore_values, threshold, na.rm = TRUE)
  return(ifelse(zscore_values > threshold_value, "High", "Low"))
}

# Method 2: ssGSEA approach (simplified version since GSVA kept failing)
calculate_ssgsea_group <- function(df, genes, threshold = 0.50) {
  print(paste("Using simplified ssGSEA for", paste(genes, collapse="+")))
  # Just using mean z-score as a fallback - Daniel said this was OK for comparison
  combined_score <- rowMeans(scale(df[genes]), na.rm = TRUE)
  threshold_value <- quantile(combined_score, threshold, na.rm = TRUE)
  return(ifelse(combined_score > threshold_value, "High", "Low"))
}

# Survival analysis using base R - ggsurvplot was too complicated
run_survival_analysis_baseR <- function(df, group_col, method_name, signature_name, endpoint, save_folder) {
  
  print(paste("    Running", method_name, "for", endpoint))
  
  # Set up column names
  time_col <- paste0(endpoint, ".time")
  event_col <- endpoint
  
  # Check I have the right columns
  required_cols <- c(time_col, event_col, group_col)
  if(!all(required_cols %in% colnames(df))) {
    print(paste("    Missing columns:", paste(required_cols[!required_cols %in% colnames(df)], collapse=", ")))
    return(NULL)
  }
  
  # Clean up the data
  clean_data <- data.frame(
    time = as.numeric(df[[time_col]]),
    event = as.numeric(df[[event_col]]),
    group = as.character(df[[group_col]]),
    stringsAsFactors = FALSE
  )
  
  # Remove problematic rows
  clean_data <- clean_data[complete.cases(clean_data), ]
  clean_data <- clean_data[clean_data$time > 0, ]
  clean_data <- clean_data[clean_data$event %in% c(0, 1), ]
  clean_data <- clean_data[clean_data$group %in% c("High", "Low"), ]
  
  if(nrow(clean_data) < 10) {
    print(paste("    Not enough data:", nrow(clean_data), "rows"))
    return(NULL)
  }
  
  # Check group sizes 
  group_table <- table(clean_data$group)
  n_high <- as.numeric(group_table["High"])
  n_low <- as.numeric(group_table["Low"])
  
  if(is.na(n_high)) n_high <- 0
  if(is.na(n_low)) n_low <- 0
  
  if(n_high < 5 || n_low < 5) {
    print(paste("    Groups too small - High:", n_high, "Low:", n_low))
    return(NULL)
  }
  
  # Make sure group is a factor
  clean_data$group <- factor(clean_data$group, levels = c("Low", "High"))
  
  # Do the survival analysis
  tryCatch({
    
    # Fit survival curves for each group
    surv_high <- survfit(Surv(time, event) ~ 1, data = clean_data[clean_data$group == "High", ])
    surv_low <- survfit(Surv(time, event) ~ 1, data = clean_data[clean_data$group == "Low", ])
    
    # Log-rank test
    test_result <- survdiff(Surv(time, event) ~ group, data = clean_data)
    pval <- 1 - pchisq(test_result$chisq, df = 1)
    
    # Make filename - had to clean up special characters
    safe_sig <- gsub("[^A-Za-z0-9_]", "_", signature_name)
    safe_method <- gsub("[^A-Za-z0-9_]", "_", method_name)
    filename <- paste0(safe_sig, "_", safe_method, "_", endpoint, ".png")
    filepath <- file.path(save_folder, filename)
    
    # Create the plot using base R
    png(filepath, width = 800, height = 600, res = 150)
    
    # Plot the curves
    plot(surv_low, col = "blue", lwd = 2, 
         xlab = "Time (days)", ylab = "Survival Probability",
         main = paste(signature_name, "-", method_name, "\n", endpoint, "| p =", round(pval, 4)))
    
    lines(surv_high, col = "red", lwd = 2)
    
    # Add legend
    legend("topright", 
           legend = c(paste("Low (n =", n_low, ")"), paste("High (n =", n_high, ")")),
           col = c("blue", "red"), 
           lwd = 2)
    
    # Add p-value
    text(x = max(clean_data$time) * 0.02, y = 0.1, 
         labels = paste("Log-rank p =", round(pval, 4)), 
         adj = 0, cex = 0.9)
    
    dev.off()
    
    print(paste("    Saved:", filename))
    
    return(list(
      signature = signature_name,
      method = method_name,
      endpoint = endpoint,
      n_high = n_high,
      n_low = n_low,
      pval = pval,
      significant = pval < 0.05,
      filename = filename
    ))
    
  }, error = function(e) {
    print(paste("    Error:", e$message))
    return(NULL)
  })
}

# Main analysis loop
results <- data.frame()

for(sig_name in names(gene_signatures)) {
  genes <- gene_signatures[[sig_name]]
  
  # Check I have all the genes
  if(!all(genes %in% colnames(DF))) {
    missing <- genes[!genes %in% colnames(DF)]
    print(paste("Skipping", sig_name, "- missing genes:", paste(missing, collapse=", ")))
    next
  }
  
  print(paste("\nAnalyzing:", sig_name))
  print(paste("Genes:", paste(genes, collapse = ", ")))
  
  # Get groups using both methods
  zscore_groups <- get_zscore_group(DF, sig_name, THRESHOLD)
  if(is.null(zscore_groups)) {
    print(paste("Failed to create z-score groups for", sig_name))
    next
  }
  
  ssgsea_groups <- calculate_ssgsea_group(DF, genes, THRESHOLD)
  
  # Add to dataframe
  z_col <- paste0("task3_z_", sig_name)
  s_col <- paste0("task3_s_", sig_name)
  DF[[z_col]] <- zscore_groups
  DF[[s_col]] <- ssgsea_groups
  
  # Check group sizes
  print(paste("Z-score groups: High =", sum(zscore_groups == "High", na.rm=T), 
              "| Low =", sum(zscore_groups == "Low", na.rm=T)))
  print(paste("ssGSEA groups: High =", sum(ssgsea_groups == "High", na.rm=T), 
              "| Low =", sum(ssgsea_groups == "Low", na.rm=T)))
  
  # Test each survival endpoint
  for(endpoint in c("OS", "DSS", "PFI")) {
    print(paste("  Testing", endpoint, ":"))
    
    # Run both methods
    z_result <- run_survival_analysis_baseR(DF, z_col, "Z-score", sig_name, endpoint, zscore_folder)
    s_result <- run_survival_analysis_baseR(DF, s_col, "ssGSEA", sig_name, endpoint, ssgsea_folder)
    
    # Compare results if both worked
    if(!is.null(z_result) && !is.null(s_result)) {
      
      # Save comparison
      comp_row <- data.frame(
        Signature = sig_name,
        Genes = paste(genes, collapse = ", "),
        Endpoint = endpoint,
        ZScore_P = round(z_result$pval, 6),
        ssGSEA_P = round(s_result$pval, 6),
        ZScore_Sig = z_result$significant,
        ssGSEA_Sig = s_result$significant,
        Better_Method = ifelse(z_result$pval < s_result$pval, "Z-score", "ssGSEA"),
        P_Difference = round(z_result$pval - s_result$pval, 6),
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, comp_row)
      
      # Show comparison
      print(paste("    Z-score: p =", round(z_result$pval, 4), 
                  ifelse(z_result$significant, "(significant)", "")))
      print(paste("    ssGSEA:  p =", round(s_result$pval, 4), 
                  ifelse(s_result$significant, "(significant)", "")))
      print(paste("    Winner: ", comp_row$Better_Method))
      
    } else {
      print("    One or both methods failed")
    }
  }
}

# Save and summarize everything
if(nrow(results) > 0) {
  
  print("\n==================================================")
  print("COMPARISON RESULTS")
  print("==================================================")
  
  # Show main results
  display_results <- results[, c("Signature", "Endpoint", "ZScore_P", "ssGSEA_P", "Better_Method")]
  print(display_results)
  
  # Save detailed results
  results_file <- file.path(comparison_folder, "task3_method_comparison.csv")
  write.csv(results, results_file, row.names = FALSE)
  print(paste("\nDetailed results saved:", results_file))
  
  # Summary stats
  print("\nSUMMARY:")
  total <- nrow(results)
  z_wins <- sum(results$Better_Method == "Z-score")
  s_wins <- sum(results$Better_Method == "ssGSEA")
  both_sig <- sum(results$ZScore_Sig & results$ssGSEA_Sig)
  either_sig <- sum(results$ZScore_Sig | results$ssGSEA_Sig)
  
  print(paste("Total comparisons:", total))
  print(paste("Z-score better:", z_wins, paste0("(", round(100*z_wins/total, 1), "%)")))
  print(paste("ssGSEA better:", s_wins, paste0("(", round(100*s_wins/total, 1), "%)")))
  print(paste("Both significant:", both_sig, paste0("(", round(100*both_sig/total, 1), "%)")))
  print(paste("Either significant:", either_sig, paste0("(", round(100*either_sig/total, 1), "%)")))
  
  # Which signatures worked best
  print("\nSIGNATURE PERFORMANCE:")
  sig_summary <- results %>%
    group_by(Signature) %>%
    summarise(
      Comparisons = n(),
      Z_Significant = sum(ZScore_Sig),
      S_Significant = sum(ssGSEA_Sig),
      Total_Significant = sum(ZScore_Sig | ssGSEA_Sig),
      Success_Rate = round(100 * Total_Significant / Comparisons, 1),
      .groups = 'drop'
    ) %>%
    arrange(desc(Success_Rate))
  
  print(sig_summary)
  
  # Create a summary report file
  report_file <- file.path(comparison_folder, "task3_summary_report.txt")
  sink(report_file)
  cat("METHOD COMPARISON SUMMARY\n")
  cat("========================\n\n")
  cat("Date:", as.character(Sys.Date()), "\n")
  cat("Total patients:", nrow(DF), "\n")
  cat("Threshold used:", THRESHOLD, "(median split)\n\n")
  cat("RESULTS:\n")
  cat("Total comparisons:", total, "\n")
  cat("Z-score performed better:", z_wins, paste0("(", round(100*z_wins/total, 1), "%)"), "\n")
  cat("ssGSEA performed better:", s_wins, paste0("(", round(100*s_wins/total, 1), "%)"), "\n")
  cat("Both methods significant:", both_sig, paste0("(", round(100*both_sig/total, 1), "%)"), "\n")
  cat("Either method significant:", either_sig, paste0("(", round(100*either_sig/total, 1), "%)"), "\n\n")
  cat("SIGNATURE PERFORMANCE:\n")
  print(sig_summary)
  sink()
  
} else {
  print("No successful comparisons - something went wrong")
}

# CORRELATION ANALYSIS - Added as requested
print("\n==================================================")
print("CORRELATION ANALYSIS")
print("==================================================")

# Create folder for correlation plots
correlation_folder <- file.path(main_folder, "correlation_plots")
if (!dir.exists(correlation_folder)) {
  dir.create(correlation_folder, recursive = TRUE)
  print(paste("Created folder:", correlation_folder))
}

# Store correlation results
correlation_results <- data.frame()

for(sig_name in names(gene_signatures)) {
  genes <- gene_signatures[[sig_name]]
  
  # Check I have all the genes (same check as before)
  if(!all(genes %in% colnames(DF))) {
    missing <- genes[!genes %in% colnames(DF)]
    print(paste("Skipping correlation for", sig_name, "- missing genes:", paste(missing, collapse=", ")))
    next
  }
  
  print(paste("Creating correlation plot for:", sig_name))
  
  # Get z-score values (continuous)
  zscore_col <- paste0("zsign_", sig_name)
  if(!zscore_col %in% colnames(DF)) {
    print(paste("Can't find z-score column:", zscore_col))
    next
  }
  zscore_values <- DF[[zscore_col]]
  
  # Calculate ssGSEA values (continuous) - same method as before
  ssgsea_values <- rowMeans(scale(DF[genes]), na.rm = TRUE)
  
  # Remove rows with missing values
  complete_data <- data.frame(
    zscore = zscore_values,
    ssgsea = ssgsea_values
  )
  complete_data <- complete_data[complete.cases(complete_data), ]
  
  if(nrow(complete_data) < 10) {
    print(paste("Not enough complete data for", sig_name, ":", nrow(complete_data), "rows"))
    next
  }
  
  # Calculate correlation
  correlation <- cor(complete_data$zscore, complete_data$ssgsea, method = "pearson")
  correlation_spearman <- cor(complete_data$zscore, complete_data$ssgsea, method = "spearman")
  
  # Statistical test
  cor_test <- cor.test(complete_data$zscore, complete_data$ssgsea, method = "pearson")
  
  print(paste("  Pearson correlation:", round(correlation, 3)))
  print(paste("  Spearman correlation:", round(correlation_spearman, 3)))
  print(paste("  P-value:", round(cor_test$p.value, 6)))
  
  # Create correlation plot
  safe_sig <- gsub("[^A-Za-z0-9_]", "_", sig_name)
  filename <- paste0("correlation_", safe_sig, ".png")
  filepath <- file.path(correlation_folder, filename)
  
  png(filepath, width = 800, height = 600, res = 150)
  
  # Create the scatter plot
  plot(complete_data$zscore, complete_data$ssgsea,
       xlab = "Z-score Signature Score",
       ylab = "ssGSEA Signature Score (simplified)",
       main = paste("Correlation:", sig_name, "\nGenes:", paste(genes, collapse = ", ")),
       pch = 16, col = rgb(0, 0, 1, 0.6), cex = 0.8)
  
  # Add regression line
  abline(lm(ssgsea ~ zscore, data = complete_data), col = "red", lwd = 2)
  
  # Add correlation info
  legend("topleft", 
         legend = c(
           paste("n =", nrow(complete_data)),
           paste("Pearson r =", round(correlation, 3)),
           paste("Spearman ρ =", round(correlation_spearman, 3)),
           paste("p =", round(cor_test$p.value, 4))
         ),
         bty = "n", cex = 0.9)
  
  # Add grid for better readability
  grid(col = "lightgray", lty = "dotted")
  
  dev.off()
  
  print(paste("  Saved correlation plot:", filename))
  
  # Store results
  corr_row <- data.frame(
    Signature = sig_name,
    Genes = paste(genes, collapse = ", "),
    N_Samples = nrow(complete_data),
    Pearson_R = round(correlation, 4),
    Spearman_Rho = round(correlation_spearman, 4),
    P_Value = round(cor_test$p.value, 6),
    Significant = cor_test$p.value < 0.05,
    R_Squared = round(correlation^2, 4),
    stringsAsFactors = FALSE
  )
  
  correlation_results <- rbind(correlation_results, corr_row)
}

# Save correlation results
if(nrow(correlation_results) > 0) {
  
  print("\nCORRELATION SUMMARY:")
  print(correlation_results[, c("Signature", "N_Samples", "Pearson_R", "Spearman_Rho", "P_Value")])
  
  # Save detailed correlation results
  corr_file <- file.path(correlation_folder, "correlation_results.csv")
  write.csv(correlation_results, corr_file, row.names = FALSE)
  print(paste("\nCorrelation results saved:", corr_file))
  
  # Summary statistics
  avg_pearson <- mean(correlation_results$Pearson_R, na.rm = TRUE)
  avg_spearman <- mean(correlation_results$Spearman_Rho, na.rm = TRUE)
  sig_correlations <- sum(correlation_results$Significant, na.rm = TRUE)
  
  print(paste("\nCorrelation Summary:"))
  print(paste("Average Pearson correlation:", round(avg_pearson, 3)))
  print(paste("Average Spearman correlation:", round(avg_spearman, 3)))
  print(paste("Significant correlations:", sig_correlations, "out of", nrow(correlation_results)))
  
  # Update the summary report with correlation info
  report_file <- file.path(comparison_folder, "task3_summary_report.txt")
  sink(report_file, append = TRUE)
  cat("\n\nCORRELATION ANALYSIS:\n")
  cat("====================\n")
  cat("Average Pearson correlation:", round(avg_pearson, 3), "\n")
  cat("Average Spearman correlation:", round(avg_spearman, 3), "\n")
  cat("Significant correlations:", sig_correlations, "out of", nrow(correlation_results), "\n\n")
  cat("DETAILED CORRELATION RESULTS:\n")
  print(correlation_results)
  sink()
  
} else {
  print("No correlation plots generated - check data")
}

# Final summary
print(paste("\nMethod comparison complete!"))
print(paste("Results saved in:", main_folder))

# Count files
z_files <- list.files(zscore_folder, pattern = "\\.png$")
s_files <- list.files(ssgsea_folder, pattern = "\\.png$")
corr_files <- list.files(correlation_folder, pattern = "\\.png$")
print(paste("Generated", length(z_files), "z-score plots,", length(s_files), "ssGSEA plots, and", length(corr_files), "correlation plots"))

if(length(z_files) > 0 || length(s_files) > 0) {
  print("\nTask 3 complete:")
  print("- Method comparison done")
  print("- Plots saved using base R")
  print("- CSV file has all the p-values")
  print("- Summary report written")
  print("- Correlation analysis included")
} else {
  print("\nNo plots generated - check data")
}