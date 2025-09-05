### 02_clean_merge_COAD_primary_4genes_FINAL.R
### Clean merge for CES1 + 3 gene signature (COAD primary only)
### Includes duplicate removal + match-fix
### Author: Alicia — Internship Week 4
### -----------------------------------------------------------------------------

library(dplyr)

# --- Set working directory ---
setwd("C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA")
MAIN_DIR <- "."

# --- Load datasets ---
expr <- read.csv(file.path(MAIN_DIR, "TOIL_RSEM_TPM.COAD.csv"), row.names = 1, check.names = FALSE)
colnames(expr) <- gsub("_", "-", colnames(expr))  # TCGA format
clin <- read.csv(file.path(MAIN_DIR, "TOIL_clinical.COAD.csv"))
sample_type <- read.table(file.path(MAIN_DIR, "TCGA_ALL_sample_type.tsv"), sep = "\t", header = TRUE)

# --- Remove ENSEMBL version numbers ---
rownames(expr) <- sub("\\..*", "", rownames(expr))
colnames(expr) <- substr(colnames(expr), 1, 15)  # Trim sample IDs to TCGA-XX-YYYY-XX

# --- Filter to Primary Tumor samples only ---
primary_ids <- sample_type %>%
  filter(sample_type == "Primary Tumor") %>%
  pull(sample) %>%
  substr(1, 15)

# --- Prefer '01' samples if duplicates exist ---
primary_ids <- primary_ids[substr(primary_ids, 14, 15) == "01"]
primary_ids <- unique(primary_ids)

# --- Subset expression matrix ---
expr_primary <- expr[, intersect(colnames(expr), primary_ids)]

# --- Extract genes of interest ---
genes_of_interest <- c(
  "ENSG00000198848",  # CES1
  "ENSG00000110090",  # CPT1A
  "ENSG00000092051",  # ACSL3
  "ENSG00000108091"   # MGLL
)
expr_subset <- expr_primary[rownames(expr_primary) %in% genes_of_interest, ]

# --- Transpose + rename ---
expr_transposed <- as.data.frame(t(expr_subset))
colnames(expr_transposed) <- c("CES1", "CPT1A", "ACSL3", "MGLL")
expr_transposed$sampleID <- rownames(expr_transposed)

# --- Prepare clinical data ---
clin$sampleID <- substr(clin$sample, 1, 15)  # Trim to TCGA ID format
clin <- clin[!duplicated(clin$sampleID), ]   # Remove clinical duplicates
clin_primary <- clin %>%
  filter(histological_type == "Colon Adenocarcinoma") %>%
  filter(sampleID %in% expr_transposed$sampleID)

# --- Final matched expression ---
expr_matched <- expr_transposed[expr_transposed$sampleID %in% clin_primary$sampleID, ]

# --- Merge ---
merged_data <- merge(clin_primary, expr_matched, by = "sampleID")

# --- Final checks ---
cat("✅ Unique primary samples:", length(unique(expr_matched$sampleID)), "\n")
cat("✅ Clinical rows after filter:", nrow(clin_primary), "\n")
cat("✅ Merged dataset rows:", nrow(merged_data), "\n")

# --- Save ---
output_csv <- "merged_COAD_primary_4genes_vFINAL.csv"
output_rdata <- "merged_COAD_primary_4genes_vFINAL.RData"

write.csv(merged_data, file = file.path(MAIN_DIR, output_csv), row.names = FALSE)
save(merged_data, file = file.path(MAIN_DIR, output_rdata))

cat("✅ Cleaned merged dataset saved as:\n - ", output_csv, "\n - ", output_rdata, "\n")
