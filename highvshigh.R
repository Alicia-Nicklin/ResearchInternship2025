# High-High analysis - Task 4 from the internship


library(survival)
library(dplyr)

# Load data as usual
setwd("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA")
DF <- read.csv("merged_COAD_primary_4genes_vFINAL.csv")

# Make folders - getting better at organizing everything
main_folder <- "KM_task4_high_high_analysis"
zscore_folder <- file.path(main_folder, "zscore_plots")
ssgsea_folder <- file.path(main_folder, "ssgsea_plots") 
results_folder <- file.path(main_folder, "results")

for(folder in c(main_folder, zscore_folder, ssgsea_folder, results_folder)) {
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
    cat("Created folder:", folder, "\n")
  }
}

# Settings - using CES1 and CPT1A like we discussed
GENES <- c("CES1", "CPT1A")
THRESHOLD <- 0.75  # 75th percentile - more selective than median

cat("High-High vs Others analysis\n")
cat("Using:", paste(GENES, collapse = " + "), "\n")
cat("Threshold:", THRESHOLD, "(75th percentile - only the highest expressers)\n")
cat("Data:", nrow(DF), "patients\n\n")

# Method 1: Z-score approach
cat("Creating z-score High-High groups...\n")

# Calculate z-scores for each gene
ces1_z <- scale(as.numeric(DF$CES1))[,1]
cpt1a_z <- scale(as.numeric(DF$CPT1A))[,1]

# Set thresholds (75th percentile for each gene)
ces1_thresh <- quantile(ces1_z, THRESHOLD, na.rm = TRUE)
cpt1a_thresh <- quantile(cpt1a_z, THRESHOLD, na.rm = TRUE)

# Create High-High group - patients must be high in BOTH genes
DF$zscore_group <- ifelse((ces1_z > ces1_thresh) & (cpt1a_z > cpt1a_thresh), 
                          "High-High", "Others")

# Method 2: ssGSEA approach (simplified version since real ssGSEA was problematic)
cat("Creating ssGSEA High-High groups...\n")

# Create different types of metabolic scores
metabolic_score <- (ces1_z + cpt1a_z) / 2  # Simple average
drug_metab_score <- ces1_z * 0.7 + cpt1a_z * 0.3  # Weighted towards CES1

# Set thresholds for both scores
metab_thresh <- quantile(metabolic_score, THRESHOLD, na.rm = TRUE)
drug_thresh <- quantile(drug_metab_score, THRESHOLD, na.rm = TRUE)

# High-High group needs high scores in BOTH pathways
DF$ssgsea_group <- ifelse((metabolic_score > metab_thresh) & (drug_metab_score > drug_thresh), 
                          "High-High", "Others")

# Check group sizes
z_high <- sum(DF$zscore_group == "High-High", na.rm = TRUE)
z_others <- sum(DF$zscore_group == "Others", na.rm = TRUE)
s_high <- sum(DF$ssgsea_group == "High-High", na.rm = TRUE)
s_others <- sum(DF$ssgsea_group == "Others", na.rm = TRUE)

cat("Z-score groups: High-High =", z_high, "| Others =", z_others, "\n")
cat("ssGSEA groups: High-High =", s_high, "| Others =", s_others, "\n\n")

# Survival analysis function - similar to what I used before
run_high_high_analysis <- function(df, group_col, method_name, endpoint, save_folder) {
  
  cat("  Analyzing", method_name, "-", endpoint, "...")
  
  # Set up columns
  time_col <- paste0(endpoint, ".time")
  event_col <- endpoint
  
  # Check I have the right columns
  if(!all(c(time_col, event_col, group_col) %in% colnames(df))) {
    cat(" FAILED (missing columns)\n")
    return(NULL)
  }
  
  # Clean up data
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
  clean_data <- clean_data[clean_data$group %in% c("High-High", "Others"), ]
  
  if(nrow(clean_data) < 10) {
    cat(" FAILED (not enough data)\n")
    return(NULL)
  }
  
  # Check group sizes
  group_table <- table(clean_data$group)
  n_high <- as.numeric(group_table["High-High"])
  n_others <- as.numeric(group_table["Others"])
  
  if(is.na(n_high)) n_high <- 0
  if(is.na(n_others)) n_others <- 0
  
  if(n_high < 5 || n_others < 5) {
    cat(" FAILED (groups too small)\n")
    return(NULL)
  }
  
  # Make group a factor
  clean_data$group <- factor(clean_data$group, levels = c("Others", "High-High"))
  
  # Do the survival analysis
  tryCatch({
    
    # Fit curves for each group
    surv_high <- survfit(Surv(time, event) ~ 1, data = clean_data[clean_data$group == "High-High", ])
    surv_others <- survfit(Surv(time, event) ~ 1, data = clean_data[clean_data$group == "Others", ])
    
    # Log-rank test
    test_result <- survdiff(Surv(time, event) ~ group, data = clean_data)
    pval <- 1 - pchisq(test_result$chisq, df = 1)
    
    # Make filename
    safe_method <- gsub("[^A-Za-z0-9_]", "_", method_name)
    filename <- paste0("CES1_CPT1A_high_high_", safe_method, "_", endpoint, ".png")
    filepath <- file.path(save_folder, filename)
    
    # Create plot
    png(filepath, width = 800, height = 600, res = 150)
    
    # Plot the curves
    plot(surv_others, col = "blue", lwd = 2, 
         xlab = "Time (days)", ylab = "Survival Probability",
         main = paste("CES1 + CPT1A High-High vs Others -", method_name, 
                      "\n", endpoint, "| p =", round(pval, 4)))
    
    lines(surv_high, col = "red", lwd = 2)
    
    # Add legend
    legend("topright", 
           legend = c(paste("Others (n =", n_others, ")"), 
                      paste("High-High (n =", n_high, ")")),
           col = c("blue", "red"), 
           lwd = 2)
    
    # Add p-value
    text(x = max(clean_data$time) * 0.02, y = 0.1, 
         labels = paste("Log-rank p =", round(pval, 4)), 
         adj = 0, cex = 0.9)
    
    dev.off()
    
    cat(" Done - p =", round(pval, 4), "\n")
    
    return(list(
      method = method_name,
      endpoint = endpoint,
      n_high = n_high,
      n_others = n_others,
      pval = pval,
      significant = pval < 0.05,
      filename = filename
    ))
    
  }, error = function(e) {
    cat(" FAILED (", e$message, ")\n")
    return(NULL)
  })
}

# Run the analyses
results <- data.frame()

cat("Running High-High survival analyses:\n")

# Test each endpoint for both methods
for(endpoint in c("OS", "DSS", "PFI")) {
  
  # Z-score method
  z_result <- run_high_high_analysis(DF, "zscore_group", "Z-score", endpoint, zscore_folder)
  
  # ssGSEA method
  s_result <- run_high_high_analysis(DF, "ssgsea_group", "ssGSEA", endpoint, ssgsea_folder)
  
  # Save results
  if(!is.null(z_result)) {
    results <- rbind(results, data.frame(
      Method = z_result$method,
      Endpoint = z_result$endpoint,
      High_High_N = z_result$n_high,
      Others_N = z_result$n_others,
      P_Value = round(z_result$pval, 6),
      Significant = z_result$significant,
      Plot_File = z_result$filename,
      stringsAsFactors = FALSE
    ))
  }
  
  if(!is.null(s_result)) {
    results <- rbind(results, data.frame(
      Method = s_result$method,
      Endpoint = s_result$endpoint,
      High_High_N = s_result$n_high,
      Others_N = s_result$n_others,
      P_Value = round(s_result$pval, 6),
      Significant = s_result$significant,
      Plot_File = s_result$filename,
      stringsAsFactors = FALSE
    ))
  }
}

# Show and save results
if(nrow(results) > 0) {
  
  cat("\n==================================================\n")
  cat("HIGH-HIGH ANALYSIS RESULTS\n")
  cat("==================================================\n\n")
  
  # Show main results
  print(results[, c("Method", "Endpoint", "High_High_N", "P_Value", "Significant")])
  
  # Save detailed results
  results_file <- file.path(results_folder, "task4_high_high_results.csv")
  write.csv(results, results_file, row.names = FALSE)
  cat("\nResults saved:", results_file, "\n")
  
  # Summary
  cat("\nSUMMARY:\n")
  total <- nrow(results)
  significant <- sum(results$Significant)
  z_sig <- sum(results$Method == "Z-score" & results$Significant)
  s_sig <- sum(results$Method == "ssGSEA" & results$Significant)
  
  cat("Total analyses:", total, "\n")
  cat("Significant results:", significant, paste0("(", round(100*significant/total, 1), "%)"), "\n")
  cat("Z-score significant:", z_sig, "\n")
  cat("ssGSEA significant:", s_sig, "\n")
  
  # Best result
  if(significant > 0) {
    best_idx <- which.min(results$P_Value)
    best <- results[best_idx, ]
    cat("\nBest result:", best$Method, "-", best$Endpoint, 
        "| p =", best$P_Value, "\n")
  }
  
  # Write summary report
  report_file <- file.path(results_folder, "task4_summary_report.txt")
  sink(report_file)
  cat("HIGH-HIGH vs OTHERS ANALYSIS SUMMARY\n")
  cat("====================================\n\n")
  cat("Gene pair: CES1 + CPT1A\n")
  cat("Threshold: 75th percentile\n")
  cat("Approach: High-High (both genes high) vs Others\n")
  cat("Date:", as.character(Sys.Date()), "\n\n")
  cat("RESULTS:\n")
  print(results)
  cat("\nSUMMARY:\n")
  cat("Total analyses:", total, "\n")
  cat("Significant results:", significant, "\n")
  sink()
  
} else {
  cat("\nNo successful analyses - something went wrong\n")
}

# Final summary
cat("\nHigh-High analysis complete!\n")
cat("Results saved in:", main_folder, "\n")

# Count files
z_files <- list.files(zscore_folder, pattern = "\\.png$")
s_files <- list.files(ssgsea_folder, pattern = "\\.png$")
cat("Generated", length(z_files), "z-score plots and", length(s_files), "ssGSEA plots\n")

if(length(z_files) > 0 || length(s_files) > 0) {
  cat("\nTask 4 complete:\n")
  cat("- High-High vs Others analysis done\n")
  cat("- Both methods tested\n") 
  cat("- Plots saved using base R\n")
  cat("- CSV and summary files created\n")
}

# Bonus: count all my survival plots from the whole project
# Set main project folder
MAIN_DIR <- "C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA"

# Find all plot files
all_plots <- list.files(MAIN_DIR, pattern = "\\.(tiff|tif|png)$", recursive = TRUE, full.names = TRUE)

# Filter for KM survival plots
km_plots <- all_plots[grepl("KM_", basename(all_plots))]

# Show total
cat("Total survival plots from entire project:", length(km_plots), "\n")
if(length(km_plots) > 0) {
  cat("First few:\n")
  head(km_plots, 5)
}