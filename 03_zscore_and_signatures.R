### TASK 3 — Normalize Gene Expression as Z-scores
### Author: Alicia
### Internship Project – Week 2
### Input: merged_COAD_primary_4genes_vFINAL.csv
### Output: merged_COAD_primary_4genes_vFINAL.csv (overwritten with added columns)
### ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
  library(survminer)
})

setwd("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA")

# Define correct input/output filenames (final version)
input_file <- "merged_COAD_primary_4genes_vFINAL.csv"
output_file <- "merged_COAD_primary_4genes_vFINAL.csv"

# Load final merged dataset
DF <- read.csv(input_file, stringsAsFactors = FALSE)

# Check required genes
required_genes <- c("CES1", "CPT1A", "ACSL3", "MGLL")
missing_genes <- setdiff(required_genes, names(DF))
if (length(missing_genes) > 0) stop("❌ Missing genes: ", paste(missing_genes, collapse=", "))

# Convert gene columns to numeric
DF[required_genes] <- lapply(DF[required_genes], function(x) as.numeric(as.character(x)))

# Compute individual z-scores
for (gene in required_genes) {
  z_col <- paste0("z_", gene)
  DF[[z_col]] <- scale(DF[[gene]], center = TRUE, scale = TRUE)[, 1]
  message("✅ Z-scored: ", gene, " → ", z_col)
}

# Define all signature combinations
signature_list <- list(
  CES1_CPT1A = c("CES1", "CPT1A"),
  CES1_MGLL_ACSL3 = c("CES1", "MGLL", "ACSL3"),
  CES1_MGLL_CPT1A = c("CES1", "MGLL", "CPT1A"),
  CES1_MGLL_CPT1A_ACSL3 = c("CES1", "MGLL", "CPT1A", "ACSL3"),
  CES1_CPT1A_ACSL3 = c("CES1", "CPT1A", "ACSL3")
)

# Compute z-score signatures + stratifications
for (sig_name in names(signature_list)) {
  genes <- signature_list[[sig_name]]
  z_data <- DF[, genes, drop = FALSE]
  z_data <- z_data[complete.cases(z_data), ]
  DF <- DF[complete.cases(DF[, genes]), ]
  
  z_mat <- scale(z_data)
  z_avg <- rowMeans(z_mat)
  
  sig_col <- paste0("zsign_", sig_name)
  DF[[sig_col]] <- z_avg
  
  qs <- quantile(z_avg, probs = c(0.25, 0.33, 0.50, 0.66, 0.75, 0.90), na.rm = TRUE)
  
  DF[[paste0(sig_col, "_25p")]] <- ifelse(z_avg > qs[1], "HIGH", "LOW")
  DF[[paste0(sig_col, "_50p")]] <- ifelse(z_avg > qs[3], "HIGH", "LOW")
  DF[[paste0(sig_col, "_75p")]] <- ifelse(z_avg > qs[5], "HIGH", "LOW")
  DF[[paste0(sig_col, "_66p")]] <- ifelse(z_avg > qs[4], "HIGH", "LOW")
  DF[[paste0(sig_col, "_90p")]] <- ifelse(z_avg > qs[6], "HIGH", "LOW")
  
  DF[[paste0(sig_col, "_25p_75p")]] <- NA
  DF[[paste0(sig_col, "_25p_75p")]][z_avg < qs[1]] <- "LOW"
  DF[[paste0(sig_col, "_25p_75p")]][z_avg > qs[5]] <- "HIGH"
  
  DF[[paste0(sig_col, "_33p_66p")]] <- NA
  DF[[paste0(sig_col, "_33p_66p")]][z_avg < qs[2]] <- "LOW"
  DF[[paste0(sig_col, "_33p_66p")]][z_avg > qs[4]] <- "HIGH"
}

# Save final output
write_csv(DF, output_file)
cat("\n✅ Final output saved to:", output_file, "\n")
