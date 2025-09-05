# 04b_signature_survival_plots.R — Follow-up threshold filtering + separate output

setwd("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA")
MAIN_DIR <- "."

# --- 1) Load packages ---
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# --- 2) Load dataset (raw) ---
raw_DF <- read.csv("merged_COAD_primary_4genes_vFINAL.csv", stringsAsFactors = FALSE)

# --- 3) Define signature combinations ---
signature_list <- list(
  CES1_CPT1A            = c("CES1", "CPT1A"),
  CES1_MGLL_ACSL3       = c("CES1", "MGLL", "ACSL3"),
  CES1_MGLL_CPT1A       = c("CES1", "MGLL", "CPT1A"),
  CES1_MGLL_CPT1A_ACSL3 = c("CES1", "MGLL", "CPT1A", "ACSL3"),
  CES1_CPT1A_ACSL3      = c("CES1", "CPT1A", "ACSL3")
)

# --- 4) Define thresholds ---
simple_q       <- c(0.25, 0.50, 0.75)
central_ranges <- list(c(0.20, 0.80), c(0.40, 0.60))
extreme_ranges <- list(c(0.00, 0.33), c(0.66, 1.00))

# --- 5) Run for each follow-up threshold ---
for (fup in c(30, 90)) {
  
  DF <- raw_DF %>% filter(OS.time >= fup, DSS.time >= fup, PFI.time >= fup)
  message("🔹 Running KM plots for follow-up ≥ ", fup, " days — n = ", nrow(DF))
  
  # --- Compute signatures + stratification ---
  # --- Use precomputed signatures (already in DF) ---
  for (sig in names(signature_list)) {
    sig_col <- paste0("zsign_", sig)
    if (!sig_col %in% colnames(DF)) {
      stop("❌ Missing expected column: ", sig_col)
    }
    
    for (q in simple_q) {
      col <- paste0("zsign_", sig, "_", q * 100, "p")
      cut <- quantile(DF[[paste0("zsign_", sig)]], q, na.rm = TRUE)
      DF[[col]] <- ifelse(DF[[paste0("zsign_", sig)]] > cut, "HIGH", "LOW")
    }
    
    for (rng in central_ranges) {
      col <- paste0("zsign_", sig, "_", rng[1]*100, "p_", rng[2]*100, "p")
      lo  <- quantile(DF[[paste0("zsign_", sig)]], rng[1], na.rm = TRUE)
      hi  <- quantile(DF[[paste0("zsign_", sig)]], rng[2], na.rm = TRUE)
      DF[[col]] <- NA
      DF[[col]][DF[[paste0("zsign_", sig)]] < lo] <- "LOW"
      DF[[col]][DF[[paste0("zsign_", sig)]] > hi] <- "HIGH"
    }
    
    for (rng in extreme_ranges) {
      col <- paste0("zsign_", sig, "_ext_", rng[1]*100, "p_", rng[2]*100, "p")
      lo  <- quantile(DF[[paste0("zsign_", sig)]], rng[1], na.rm = TRUE)
      hi  <- quantile(DF[[paste0("zsign_", sig)]], rng[2], na.rm = TRUE)
      DF[[col]] <- NA
      DF[[col]][DF[[paste0("zsign_", sig)]] <= lo] <- "LOW"
      DF[[col]][DF[[paste0("zsign_", sig)]] >= hi] <- "HIGH"
    }
  }
  
  # --- Set output directory & results table ---
  out_dir <- paste0("KM_signature_plots_min", fup)
  dir.create(out_dir, showWarnings = FALSE)
  results <- data.frame()
  
  # --- Generate plots ---
  for (sig in names(signature_list)) {
    base_col   <- paste0("zsign_", sig)
    strat_cols <- grep(paste0("^", base_col), names(DF), value = TRUE)
    
    for (strat in strat_cols) {
      vals <- DF[[strat]]
      if (all(is.na(vals)) || length(unique(na.omit(vals))) < 2) next
      
      tmp <- DF %>%
        filter(!is.na(OS), !is.na(OS.time), !is.na(.data[[strat]])) %>%
        mutate(Group = factor(.data[[strat]], levels = c("LOW", "HIGH")))
      
      if (nrow(tmp) == 0 || length(unique(tmp$Group)) < 2) next
      if (sum(tmp$Group == "LOW") < 2 || sum(tmp$Group == "HIGH") < 2) next
      
      fit  <- survfit(Surv(OS.time, OS) ~ Group, data = tmp)
      pval <- surv_pvalue(fit, data = tmp)$pval
      
      legend_labs <- c(
        paste0("Low  (n=", sum(tmp$Group == "LOW"), ")"),
        paste0("High (n=", sum(tmp$Group == "HIGH"), ")")
      )
      
      p <- ggsurvplot(
        fit, data = tmp,
        palette     = c("black", "red"),
        legend.labs = legend_labs,
        pval        = TRUE,
        xlab        = "Days",
        ylab        = "Overall Survival (%)",
        fun         = function(y) y * 100,
        title       = paste(
          "Signature:", paste(gsub("z_", "", signature_list[[sig]]), collapse = " + "),
          "\nStratification:", gsub("^zsign_", "", strat),
          paste0("\nFollow-up ≥ ", fup, " days")
        )
      )
      
      out_file <- file.path(out_dir, paste0("KM_", sig, "_", strat, "_min", fup, "days.png"))
      ggsave(out_file, plot = p$plot, width = 6, height = 5)
      
      results <- rbind(results, data.frame(
        Followup_Days  = fup,
        Signature      = sig,
        Stratification = strat,
        p_value        = round(pval, 4),
        N              = nrow(tmp),
        n_LOW          = sum(tmp$Group == "LOW"),
        n_HIGH         = sum(tmp$Group == "HIGH")
      ))
    }
  }
  
  # --- Save result summary ---
  write.csv(results, paste0("KM_pvalue_summary_filtered_min", fup, ".csv"), row.names = FALSE)
}
message("✅ All KM plots for all thresholds saved.")
