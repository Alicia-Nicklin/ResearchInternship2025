# 04c_check_zscore_signatures.R — Validate Z-Score Signature Calculations
# Author: Alicia
# Internship Week 2 — QA Script

# --- Setup ---
setwd("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA")

library(dplyr)
library(ggplot2)
library(reshape2)

# --- Load final dataset ---
DF <- read.csv("merged_COAD_primary_4genes_vFINAL.csv", stringsAsFactors = FALSE)

# --- Define signature sets ---
signature_list <- list(
  CES1_CPT1A            = c("CES1", "CPT1A"),
  CES1_MGLL_ACSL3       = c("CES1", "MGLL", "ACSL3"),
  CES1_MGLL_CPT1A       = c("CES1", "MGLL", "CPT1A"),
  CES1_MGLL_CPT1A_ACSL3 = c("CES1", "MGLL", "CPT1A", "ACSL3"),
  CES1_CPT1A_ACSL3      = c("CES1", "CPT1A", "ACSL3")
)

# --- 1. Manual inspection: show 5 samples per signature ---
cat("\n🔍 Inspecting 5 random samples for each signature:\n")
for (sig in names(signature_list)) {
  genes <- signature_list[[sig]]
  sig_col <- paste0("zsign_", sig)
  
  if (sig_col %in% names(DF)) {
    cat("\n---", sig_col, "---\n")
    print(DF %>%
            select(Sample, all_of(genes), all_of(sig_col)) %>%
            slice_sample(n = 5))
  } else {
    cat("❌ Signature column not found:", sig_col, "\n")
  }
}

# --- 2. Recalculate and compare signature scores ---
cat("\n✅ Verifying mean z-score calculations match stored values:\n")
for (sig in names(signature_list)) {
  genes <- signature_list[[sig]]
  sig_col <- paste0("zsign_", sig)
  
  if (all(genes %in% names(DF)) && sig_col %in% names(DF)) {
    mat <- scale(DF[, genes])
    comp <- rowMeans(mat)
    ok <- all.equal(round(comp, 6), round(DF[[sig_col]], 6))
    cat(sprintf("✔ %s: %s\n", sig_col, ifelse(isTRUE(ok), "MATCH", "MISMATCH ❌")))
  } else {
    cat(sprintf("❌ Missing genes or column for: %s\n", sig))
  }
}

# --- 3. Plot raw z-score distributions ---
cat("\n📊 Plotting histograms and boxplots of signature z-scores...\n")
z_cols <- grep("^zsign_[A-Z]", names(DF), value = TRUE)
z_long <- DF %>%
  select(Sample, all_of(z_cols)) %>%
  mutate(across(-Sample, as.numeric)) %>%
  melt(id.vars = "Sample")

# Histograms
ggplot(z_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms of Signature Z-Scores", x = "Z-score", y = "Count")

# Boxplots
ggplot(z_long, aes(x = variable, y = value)) +
  geom_boxplot(fill = "palegreen") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Boxplots of Signature Z-Scores", x = "Signature", y = "Z-score")

# --- 4. Check all signature genes are in the dataset ---
cat("\n🔎 Verifying required genes exist in DF:\n")
for (sig in names(signature_list)) {
  genes <- signature_list[[sig]]
  missing <- setdiff(genes, names(DF))
  if (length(missing) == 0) {
    cat("✔", sig, ": all genes present\n")
  } else {
    cat("❌", sig, ": missing genes →", paste(missing, collapse = ", "), "\n")
  }
}

cat("\n✅ Z-score signature QA complete.\n")
