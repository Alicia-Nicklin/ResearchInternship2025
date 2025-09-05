# --- 0) Load packages ---
library(survival)
library(survminer)
library(dplyr)

# --- 1) Load merged dataset ---
df <- read.csv("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA/merged_COAD_primary_4genes_vFINAL.csv",
               stringsAsFactors = FALSE,
               check.names = TRUE)

# --- 2) Compute composite signature from specified genes ---
compute_sig <- function(df, genes, out.col = "sig") {
  missing <- setdiff(genes, names(df))
  if (length(missing) > 0) {
    stop("❌ Missing genes in dataset: ", paste(missing, collapse = ", "))
  }
  mat <- df[, genes, drop = FALSE]
  complete_rows <- complete.cases(mat)
  mat_scaled <- scale(mat[complete_rows, , drop = FALSE])
  row_means  <- rowMeans(mat_scaled, na.rm = TRUE)
  df[[out.col]] <- NA
  df[[out.col]][complete_rows] <- row_means
  return(df)
}

# --- 3) Stratify a numeric vector into common quantile schemes ---
stratify_vec <- function(x) {
  x <- as.numeric(x)
  out <- list()
  for (p in c(0.25, 0.50, 0.75)) {
    nm  <- paste0("q", p * 100)
    thr <- quantile(x, p, na.rm = TRUE)
    out[[nm]] <- ifelse(x > thr, "HIGH", "LOW")
  }
  for (rng in list(c(0.20, 0.80), c(0.40, 0.60))) {
    nm <- paste0("c", rng[1]*100, "_", rng[2]*100)
    lo <- quantile(x, rng[1], na.rm = TRUE)
    hi <- quantile(x, rng[2], na.rm = TRUE)
    v <- rep(NA, length(x))
    v[x < lo] <- "LOW"
    v[x > hi] <- "HIGH"
    out[[nm]] <- v
  }
  for (rng in list(c(0.00, 0.33), c(0.66, 1.00))) {
    nm <- paste0("e", rng[1]*100, "_", rng[2]*100)
    lo <- quantile(x, rng[1], na.rm = TRUE)
    hi <- quantile(x, rng[2], na.rm = TRUE)
    v <- rep(NA, length(x))
    v[x <= lo] <- "LOW"
    v[x >= hi] <- "HIGH"
    out[[nm]] <- v
  }
  return(out)
}

# --- 4) Plot + save KM curve ---
plot_km <- function(df_sub, genes_used, strat_name, out_base_dir, surv_type) {
  if (!all(c("Group", paste0(surv_type, ".time"), surv_type) %in% colnames(df_sub))) return()
  
  fit <- survfit(Surv(df_sub[[paste0(surv_type, ".time")]], df_sub[[surv_type]]) ~ Group, data = df_sub)
  
  # Detect if pval should be safely calculated
  safe_for_pval <- tryCatch({
    survminer::surv_pvalue(fit, data = df_sub)
    TRUE
  }, error = function(e) {
    message("⚠️ Skipping p-value due to error: ", e$message)
    FALSE
  })
  
  legend_labs <- c(
    paste0("Low  (n=", sum(df_sub$Group == "LOW"), ")"),
    paste0("High (n=", sum(df_sub$Group == "HIGH"), ")")
  )
  
  plot_title <- paste(
    "Signature:", paste(genes_used, collapse = " + "),
    "\nStratification:", strat_name,
    paste0("\nSurvival Type: ", surv_type)
  )
  
  if (safe_for_pval) {
    p <- ggsurvplot(
      fit, data = df_sub,
      palette     = c("black", "red"),
      legend.labs = legend_labs,
      pval        = TRUE,
      xlab        = "Days",
      ylab        = paste0(surv_type, " (%)"),
      fun         = function(y) y * 100,
      title       = plot_title
    )
  } else {
    p <- ggsurvplot(
      fit, data = df_sub,
      palette     = c("black", "red"),
      legend.labs = legend_labs,
      pval        = FALSE,
      xlab        = "Days",
      ylab        = paste0(surv_type, " (%)"),
      fun         = function(y) y * 100,
      title       = plot_title
    )
  }
  
  out_dir <- file.path(out_base_dir, surv_type)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fname <- file.path(
    out_dir,
    paste0(surv_type, "_KM_", paste(genes_used, collapse = "_"), "_", strat_name, ".tiff")
  )
  ggsave(fname, plot = p$plot, width = 6, height = 5)
}

# --- 5) Signature list ---
signature_list <- list(
  CES1_CPT1A                  = c("CES1", "CPT1A"),
  CES1_CPT1A_ACSL3            = c("CES1", "CPT1A", "ACSL3"),
  CES1_CPT1A_ACSL3_MGLL       = c("CES1", "CPT1A", "ACSL3", "MGLL"),
  CES1_ACSL3_MGLL             = c("CES1", "ACSL3", "MGLL"),
  CES1_MGLL                   = c("CES1", "MGLL"),
  CPT1A_ACSL3_MGLL            = c("CPT1A", "ACSL3", "MGLL")
)

# --- 6) Output folder ---
out_base_dir <- "KM_signature_plots"

# --- 7) Main loop ---
surv_types <- c("OS", "DSS", "PFI")

for (surv_type in surv_types) {
  for (sig_name in names(signature_list)) {
    genes <- signature_list[[sig_name]]
    df2   <- compute_sig(df, genes, out.col = "sig")
    df2   <- df2[!is.na(df2$sig), ]
    
    if (nrow(df2) == 0) next
    
    strat <- stratify_vec(df2$sig)
    
    for (scheme in names(strat)) {
      group <- strat[[scheme]]
      keep <- !is.na(group)
      group <- group[keep]
      df_sub <- df2[keep, ]
      
      if (length(group) != nrow(df_sub)) {
        cat("❌ Mismatch: group=", length(group), " vs df_sub=", nrow(df_sub), "\n")
        next
      }
      
      if (length(unique(group)) < 2) next
      
      df_sub$Group <- factor(group, levels = c("LOW", "HIGH"))
      
      plot_km(df_sub,
              genes_used   = genes,
              strat_name   = scheme,
              out_base_dir = out_base_dir,
              surv_type    = surv_type)
    }
  }
}

message("✅ All survival plots saved to folders in ", out_base_dir, "/")
