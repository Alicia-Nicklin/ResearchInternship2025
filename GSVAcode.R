# ssGSEA survival analysis - finally got this working!
# Had so many issues with GSVA installation 

# Install packages 
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GSVA", "GSEABase", "survival", "survminer"), ask = FALSE)

library(GSVA)
library(GSEABase)
library(survival)
library(survminer)
library(ggplot2)

# Load my data 
print("Loading data...")

setwd("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA")
df <- read.csv("merged_COAD_primary_4genes_vFINAL.csv")

# Alternative ways to load if the above doesn't work
# df <- read.csv("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA/merged_COAD_primary_4genes_vFINAL.csv")
# df <- read.csv(file.choose())  # opens dialog

print(paste("Data loaded:", nrow(df), "rows x", ncol(df), "columns"))
print("First few column names:")
print(head(colnames(df), 20))

# Check I have the survival columns I need
required_columns <- c("sampleID", "OS.time", "OS", "DSS.time", "DSS", "PFI.time", "PFI")
missing_columns <- required_columns[!required_columns %in% colnames(df)]
if(length(missing_columns) > 0) {
  print("Missing columns:")
  print(missing_columns)
} else {
  print("Good - all survival columns found")
}

# My genes of interest
genes_of_interest <- c("CES1", "CPT1A", "ACSL3", "MGLL")
available_genes <- genes_of_interest[genes_of_interest %in% colnames(df)]

print(paste("Available genes:", paste(available_genes, collapse = ", ")))

missing_genes <- genes_of_interest[!genes_of_interest %in% colnames(df)]
if(length(missing_genes) > 0) {
  print("Missing genes:")
  print(missing_genes)
}

if(length(available_genes) == 0) {
  stop("No target genes found - check column names")
}

# Prepare expression matrix for GSVA
# This part was confusing - GSVA wants genes as rows, samples as columns
expr_df <- df[, available_genes, drop = FALSE]
expr_df <- data.frame(lapply(expr_df, function(x) as.numeric(as.character(x))))
rownames(expr_df) <- df$sampleID

# Transpose for GSVA 
expr_matrix <- t(expr_df)

# Clean up NAs - GSVA is picky about this
expr_matrix <- expr_matrix[, colSums(is.na(expr_matrix)) == 0]
expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0, ]
mode(expr_matrix) <- "numeric"

print(paste("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples"))

# Prepare gene sets - GSVA needs this as a list
print("Preparing gene sets...")

gene_set_list <- list(
  CES1_signature = available_genes
  # Could add more signatures here later
)

print("Gene sets:")
print(gene_set_list)

# Run ssGSEA - this is where I had most problems
print("Running ssGSEA...")

ssgsea_scores <- NULL
ssgsea_success <- FALSE

# Try the new GSVA version first
tryCatch({
  print("Trying new GSVA method...")
  
  ssgsea_param <- ssgseaParam(
    exprData = expr_matrix, 
    geneSets = gene_set_list,
    alpha = 0.25,
    normalize = TRUE
  )
  
  ssgsea_scores <- gsva(ssgsea_param, verbose = TRUE)
  ssgsea_success <- TRUE
  print("New method worked!")
  
}, error = function(e1) {
  print(paste("New method failed:", e1$message))
  
  # Try the old way
  tryCatch({
    print("Trying legacy method...")
    
    ssgsea_scores <- gsva(
      expr_matrix, 
      gene_set_list, 
      method = "ssgsea",
      ssgsea.norm = TRUE,  # Daniel said this was important
      verbose = TRUE
    )
    
    ssgsea_success <- TRUE
    print("Legacy method worked!")
    
  }, error = function(e2) {
    print(paste("Legacy method also failed:", e2$message))
    print("Using backup calculation...")
  })
})

# Backup method if GSVA completely fails
if (!ssgsea_success || is.null(ssgsea_scores)) {
  print("Using custom calculation since GSVA failed...")
  
  # Simple ssGSEA-like calculation
  calculate_ssgsea_scores <- function(expr_matrix, gene_set_list, alpha = 0.25) {
    
    scores_matrix <- matrix(NA, nrow = length(gene_set_list), ncol = ncol(expr_matrix))
    rownames(scores_matrix) <- names(gene_set_list)
    colnames(scores_matrix) <- colnames(expr_matrix)
    
    for (set_name in names(gene_set_list)) {
      genes_in_set <- intersect(gene_set_list[[set_name]], rownames(expr_matrix))
      
      if (length(genes_in_set) == 0) {
        warning(paste("No genes found for set:", set_name))
        next
      }
      
      print(paste("Processing:", set_name, "with", length(genes_in_set), "genes"))
      
      # Calculate enrichment for each sample
      enrichment_scores <- sapply(colnames(expr_matrix), function(sample) {
        sample_expr <- expr_matrix[genes_in_set, sample, drop = FALSE]
        
        if (length(genes_in_set) == 1) {
          return(sample_expr[1, 1])
        } else {
          # Rank-based weighted mean (trying to mimic ssGSEA)
          gene_ranks <- rank(sample_expr[, 1], ties.method = "average")
          weights <- gene_ranks^alpha
          weighted_score <- sum(sample_expr[, 1] * weights) / sum(weights)
          return(weighted_score)
        }
      })
      
      scores_matrix[set_name, ] <- enrichment_scores
    }
    
    return(scores_matrix)
  }
  
  ssgsea_scores <- calculate_ssgsea_scores(expr_matrix, gene_set_list, alpha = 0.25)
  
  # Normalize like ssgsea.norm = TRUE
  ssgsea_scores <- t(scale(t(ssgsea_scores)))
}

# Check the scores look reasonable
print("Checking ssGSEA scores...")
print(paste("Score matrix:", nrow(ssgsea_scores), "x", ncol(ssgsea_scores)))
print("First 5 scores:")
print(round(ssgsea_scores[1, 1:5], 4))

score_vector <- as.vector(ssgsea_scores[1, ])
print(paste("Score range:", round(min(score_vector, na.rm = TRUE), 4), "to", round(max(score_vector, na.rm = TRUE), 4)))
print(paste("Mean:", round(mean(score_vector, na.rm = TRUE), 4)))

# Get patient-level scores
print("Extracting patient scores...")

patient_ssgsea_scores <- as.vector(ssgsea_scores[1, ])  # First gene set
names(patient_ssgsea_scores) <- colnames(ssgsea_scores)

print(paste("Got scores for", length(patient_ssgsea_scores), "patients"))

# Match with clinical data
clinical_data <- df[match(names(patient_ssgsea_scores), df$sampleID), ]
clinical_data$ssgsea_score <- patient_ssgsea_scores[clinical_data$sampleID]

# Group patients by quantiles like Daniel suggested
print("Grouping patients...")

valid_scores <- !is.na(clinical_data$ssgsea_score) & is.finite(clinical_data$ssgsea_score)
valid_score_values <- clinical_data$ssgsea_score[valid_scores]

# Calculate quantiles
q25 <- quantile(valid_score_values, 0.25, na.rm = TRUE)
q50 <- quantile(valid_score_values, 0.50, na.rm = TRUE)
q75 <- quantile(valid_score_values, 0.75, na.rm = TRUE)

print(paste("Quantiles - 25%:", round(q25, 4), "50%:", round(q50, 4), "75%:", round(q75, 4)))

# Create groups - trying both 3-way and binary splits
clinical_data$ssgsea_group <- cut(
  clinical_data$ssgsea_score,
  breaks = c(-Inf, q25, q75, Inf),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)

clinical_data$ssgsea_binary <- ifelse(
  is.na(clinical_data$ssgsea_score), NA,
  ifelse(clinical_data$ssgsea_score > q50, "High", "Low")
)

print("Group sizes:")
print(table(clinical_data$ssgsea_group, useNA = "always"))
print(table(clinical_data$ssgsea_binary, useNA = "always"))

# Survival analysis function - this got quite long
run_ssgsea_survival_analysis <- function(clinical_data, survival_type = "OS") {
  
  print(paste("Running", survival_type, "analysis..."))
  
  # Set up survival columns
  if(survival_type == "OS") {
    time_col <- "OS.time"
    event_col <- "OS"
  } else if(survival_type == "DSS") {
    time_col <- "DSS.time" 
    event_col <- "DSS"
  } else if(survival_type == "PFI") {
    time_col <- "PFI.time"
    event_col <- "PFI"
  } else {
    stop("Invalid survival type")
  }
  
  # Check columns exist
  if(!all(c(time_col, event_col) %in% colnames(clinical_data))) {
    print(paste("Missing", survival_type, "columns"))
    return(NULL)
  }
  
  # Clean data
  surv_data <- clinical_data[
    !is.na(clinical_data[[time_col]]) & 
      !is.na(clinical_data[[event_col]]) &
      !is.na(clinical_data$ssgsea_score), 
  ]
  
  if(nrow(surv_data) < 10) {
    print(paste("Not enough samples for", survival_type))
    return(NULL)
  }
  
  print(paste("Analyzing", nrow(surv_data), "samples"))
  
  # Make standard columns for plotting
  surv_data$time <- surv_data[[time_col]]
  surv_data$event <- surv_data[[event_col]]
  
  # Fit models
  print("Fitting models...")
  
  # Continuous score
  fit_continuous <- coxph(Surv(time, event) ~ ssgsea_score, data = surv_data)
  
  # Binary groups
  fit_binary <- survfit(Surv(time, event) ~ ssgsea_binary, data = surv_data)
  
  # Three groups  
  fit_groups <- survfit(Surv(time, event) ~ ssgsea_group, data = surv_data)
  
  # Log-rank tests
  logrank_binary <- survdiff(Surv(time, event) ~ ssgsea_binary, data = surv_data)
  logrank_groups <- survdiff(Surv(time, event) ~ ssgsea_group, data = surv_data)
  
  # Results
  print("Continuous model:")
  print(summary(fit_continuous))
  
  binary_pval <- 1 - pchisq(logrank_binary$chisq, df = 1)
  groups_pval <- 1 - pchisq(logrank_groups$chisq, df = 2)
  
  print(paste("Binary groups p-value:", round(binary_pval, 4)))
  print(paste("Three groups p-value:", round(groups_pval, 4)))
  
  # Make plots
  print("Making plots...")
  
  # Binary plot
  tryCatch({
    p1 <- ggsurvplot(
      fit_binary, 
      data = surv_data,
      title = paste(survival_type, "Survival by ssGSEA Score (Binary)"),
      subtitle = "ssGSEA Analysis - CES1 Signature",
      xlab = paste(survival_type, "Time (days)"),
      ylab = paste(survival_type, "Probability"),
      pval = TRUE,
      pval.method = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.3,
      legend.title = "ssGSEA Group",
      legend.labs = c("High", "Low"),
      palette = c("red", "blue"),
      ggtheme = theme_bw(base_size = 14) + theme(
        plot.title = element_text(hjust = 0.5),
        legend.background = element_blank(),
        panel.grid = element_blank()
      )
    )
    
    filename_binary <- paste0("ssgsea_", tolower(survival_type), "_binary_survival.png")
    ggsave(filename = filename_binary, plot = p1$plot, width = 12, height = 10, dpi = 300)
    print(paste("Saved binary plot:", filename_binary))
    
  }, error = function(e) {
    print(paste("Error with binary plot:", e$message))
  })
  
  # Three groups plot
  tryCatch({
    p2 <- ggsurvplot(
      fit_groups, 
      data = surv_data,
      title = paste(survival_type, "Survival by ssGSEA Score (Three Groups)"),
      subtitle = "ssGSEA Analysis - CES1 Signature", 
      xlab = paste(survival_type, "Time (days)"),
      ylab = paste(survival_type, "Probability"),
      pval = TRUE,
      pval.method = TRUE,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.3,
      legend.title = "ssGSEA Group",
      legend.labs = c("High", "Low", "Medium"),
      palette = c("red", "blue", "green"),
      ggtheme = theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
    )
    
    filename_groups <- paste0("ssgsea_", tolower(survival_type), "_groups_survival.png")
    ggsave(filename = filename_groups, plot = p2$plot, width = 12, height = 10, dpi = 300)
    print(paste("Saved groups plot:", filename_groups))
    
  }, error = function(e) {
    print(paste("Error with groups plot:", e$message))
  })
  
  # Return everything
  return(list(
    continuous_model = fit_continuous,
    binary_fit = fit_binary,
    groups_fit = fit_groups,
    binary_pval = binary_pval,
    groups_pval = groups_pval,
    n_samples = nrow(surv_data),
    data = surv_data
  ))
}

# Run all survival analyses
survival_results <- list()

# OS
if("OS.time" %in% colnames(clinical_data) && "OS" %in% colnames(clinical_data)) {
  survival_results$OS <- run_ssgsea_survival_analysis(clinical_data, "OS")
}

# DSS  
if("DSS.time" %in% colnames(clinical_data) && "DSS" %in% colnames(clinical_data)) {
  survival_results$DSS <- run_ssgsea_survival_analysis(clinical_data, "DSS")
}

# PFI
if("PFI.time" %in% colnames(clinical_data) && "PFI" %in% colnames(clinical_data)) {
  survival_results$PFI <- run_ssgsea_survival_analysis(clinical_data, "PFI")
}

# Save results
print("Saving results...")

# Patient scores with groups
ssgsea_results_df <- data.frame(
  sampleID = names(patient_ssgsea_scores),
  ssgsea_CES1_signature = patient_ssgsea_scores,
  ssgsea_group_3levels = clinical_data$ssgsea_group[match(names(patient_ssgsea_scores), clinical_data$sampleID)],
  ssgsea_group_binary = clinical_data$ssgsea_binary[match(names(patient_ssgsea_scores), clinical_data$sampleID)],
  quantile_25th = q25,
  quantile_50th = q50,
  quantile_75th = q75
)

write.csv(ssgsea_results_df, file = "ssgsea_patient_scores_and_groups.csv", row.names = FALSE)
print("Saved patient scores to: ssgsea_patient_scores_and_groups.csv")

# Summary table
survival_summary <- data.frame(
  Analysis = c("OS_binary", "OS_groups", "DSS_binary", "DSS_groups", "PFI_binary", "PFI_groups"),
  P_value = c(
    if(!is.null(survival_results$OS)) survival_results$OS$binary_pval else NA,
    if(!is.null(survival_results$OS)) survival_results$OS$groups_pval else NA,
    if(!is.null(survival_results$DSS)) survival_results$DSS$binary_pval else NA,
    if(!is.null(survival_results$DSS)) survival_results$DSS$groups_pval else NA,
    if(!is.null(survival_results$PFI)) survival_results$PFI$binary_pval else NA,
    if(!is.null(survival_results$PFI)) survival_results$PFI$groups_pval else NA
  ),
  N_samples = c(
    if(!is.null(survival_results$OS)) survival_results$OS$n_samples else NA,
    if(!is.null(survival_results$OS)) survival_results$OS$n_samples else NA,
    if(!is.null(survival_results$DSS)) survival_results$DSS$n_samples else NA,
    if(!is.null(survival_results$DSS)) survival_results$DSS$n_samples else NA,
    if(!is.null(survival_results$PFI)) survival_results$PFI$n_samples else NA,
    if(!is.null(survival_results$PFI)) survival_results$PFI$n_samples else NA
  )
)

write.csv(survival_summary, file = "ssgsea_survival_summary.csv", row.names = FALSE)
print("Saved survival summary to: ssgsea_survival_summary.csv")

# Final check
print("ssGSEA analysis complete!")
print("Generated files:")

generated_files <- list.files(pattern = "ssgsea.*\\.(png|csv)$")
for(file in generated_files) {
  print(paste("-", file))
}

print("Done!")