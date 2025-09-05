# ResearchInternship2025
# CES1 Expression and Survival Analysis in Colorectal Cancer

This repository contains code, data, and documentation from a research project exploring the prognostic relevance of CES1 expression in colorectal cancer (COAD). The project involves survival analysis, clinical subsetting, and differential gene expression analysis using R.

## Project Overview

The analysis investigates whether CES1 expression levels are associated with patient outcomes across multiple clinical endpoints. Patients are stratified into CES1 HIGH and LOW groups, and survival analysis is performed using Overall Survival (OS), Disease-Specific Survival (DSS), and Progression-Free Interval (PFI). The workflow also includes subgroup analyses based on clinical features (e.g. tumor stage) and differential expression analysis (DEA) between expression groups.

## Repository Contents

- `data/` – Input datasets (clinical and expression data)
- `output/` – Generated Kaplan–Meier plots, DEA results, and volcano plots
- `scripts/` – R scripts for each stage of the pipeline
- `KM_plots/` – Survival curves for OS, DSS, and PFI
- `DEA_results/` – Differentially expressed genes (full, upregulated, downregulated)
- `volcano/` – Volcano plots from DEA results
- `TASKS 1-10.docx` – Full step-by-step pipeline with R code
- `README.md` – Project summary and usage guide

## Setup Instructions

Install the required R packages:

```r
install.packages(c("survminer", "survival", "reshape2", "dplyr", "ggplot2", "corrplot"))
install.packages("BiocManager")
BiocManager::install("limma")

How to Run

Clone or download this repository and open the project in RStudio.

Set the MAIN_DIR variable in your script to the full path of the project directory.

Follow the code steps outlined in TASKS 1-10.docx or the scripts folder.

Run survival analysis using Kaplan–Meier curves for OS, DSS, and PFI.

Subset by clinical features (e.g. exclude Stage IV patients) and repeat the survival analysis.

Run DEA using the limma package to compare HIGH vs LOW CES1 expression.

Generate volcano plots to visualize differentially expressed genes.

| Task    | Description                                                                 |
| ------- | --------------------------------------------------------------------------- |
| Task 1  | Install necessary R and Bioconductor packages                               |
| Task 2  | Load all required libraries                                                 |
| Task 3  | Import clinical and expression datasets                                     |
| Task 4  | Stratify patients into CES1 HIGH and LOW using quantile thresholds          |
| Task 5  | Perform Kaplan–Meier analysis for Overall Survival (OS)                     |
| Task 6  | Perform Kaplan–Meier analysis for Disease-Specific Survival (DSS)           |
| Task 7  | Perform Kaplan–Meier analysis for Progression-Free Interval (PFI)           |
| Task 8  | Subset the dataset by clinical variables (e.g. tumor stage)                 |
| Task 9  | Run differential expression analysis (DEA) between HIGH and LOW CES1 groups |
| Task 10 | Create volcano plots of differentially expressed genes                      |


Example Outputs

KM_OS_0.5_v1.tiff – Kaplan–Meier OS curve at median CES1 threshold

DEA_COAD.csv – All differentially expressed genes

UP-reg_COAD.csv, DN-reg_COAD.csv – Subsets of up- and down-regulated genes

volcano_plot.pdf – Volcano plot highlighting most significant genes

Notes

Thresholds for expression stratification (e.g. 0.25, 0.5, 0.75) can be easily adjusted to explore different quantile splits.

File naming conventions reflect analysis type, threshold used, and clinical subset (if any).

For reproducibility, all code is written using base R and widely-used packages. Paths are set using MAIN_DIR for portability.

Contact

This work was developed during a 2025 research internship in cancer bioinformatics. For questions or collaboration inquiries, please reach out via the associated GitHub profile.
