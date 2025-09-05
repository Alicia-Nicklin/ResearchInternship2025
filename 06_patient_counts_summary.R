# --- 0) Load libraries ---
library(dplyr)
library(readr)

# --- 1) Load expression (after deduplication) ---
expr <- read_csv("merged_COAD_primary_4genes_vFINAL.csv")
load('merged_COAD_primary_4genes_vFINAL.RData')

# Step 1: Expression sample count
expr_count <- nrow(expr)

# Step 2: Clinical data (survival info must be non-missing)
os_count   <- sum(!is.na(expr$OS) & !is.na(expr$OS.time))
dss_count  <- sum(!is.na(expr$DSS) & !is.na(expr$DSS.time))
pfi_count  <- sum(!is.na(expr$PFI) & !is.na(expr$PFI.time))

# Step 3 + 4: After follow-up filters
followup_30 <- expr %>%
  filter(OS.time >= 30, DSS.time >= 30, PFI.time >= 30)

followup_90 <- expr %>%
  filter(OS.time >= 90, DSS.time >= 90, PFI.time >= 90)

os_30  <- sum(!is.na(followup_30$OS) & !is.na(followup_30$OS.time))
os_90  <- sum(!is.na(followup_90$OS) & !is.na(followup_90$OS.time))

# --- 5) Summary table ---
patient_summary <- data.frame(
  Stage                             = c("Expression (after deduplication)",
                                        "With OS info",
                                        "With DSS info",
                                        "With PFI info",
                                        "After follow-up ≥ 30 days (OS)",
                                        "After follow-up ≥ 90 days (OS)"),
  Patient_Count                     = c(expr_count, os_count, dss_count, pfi_count, os_30, os_90)
)

# --- 6) Save or print ---
print(patient_summary)
write_csv(patient_summary, "patient_summary_table.csv")




