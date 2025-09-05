# --- 0) Set working directory ---
setwd("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA")

# --- 1) Load libraries ---
library(dplyr)
library(stringr)
library(readr)

# --- 2) Load raw expression and clinical data ---
# NOTE: .RData files must be loaded with `load()`, not read_csv()
load("merged_COAD_primary_4genes_vFINAL.RData")  # loads object `merged_data`
expr <- merged_data

clin_file <- "TOIL_RSEM_TPM_COAD_CES1_clinical.csv"  # update if needed
clin <- read_csv(clin_file)

# --- 3) Extract patient ID and sample type from sampleID ---
expr <- expr %>%
  mutate(
    patientID = str_sub(sampleID, 1, 12),
    sample_type = str_sub(sampleID, 14, 15)
  )

# --- 4) Keep only one sample per patient (prefer "01") ---
expr_dedup <- expr %>%
  arrange(patientID, sample_type) %>%
  group_by(patientID) %>%
  slice(1) %>%
  ungroup()

# --- 5) Merge with clinical data ---
clin <- clin %>%
  mutate(patientID = str_sub(sampleID, 1, 12))  # ensure patientID column

merged <- inner_join(clin, expr_dedup, by = "patientID")

# --- 6) Check uniqueness ---
stopifnot(!anyDuplicated(merged$patientID))
stopifnot(!anyDuplicated(merged$sampleID))

# --- 7) Save cleaned merged dataset ---
write_csv(merged, "merged_COAD_primary_4genes_v4_cleaned.csv")
message("✅ Clean merged file saved as: merged_COAD_primary_4genes_v4_cleaned.csv")
