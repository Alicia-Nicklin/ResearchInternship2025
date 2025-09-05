library("survminer")
library("survival")
library("reshape2")
library("dplyr")
library("ggplot2")
library("corrplot")

# Set working directory to local DATA folder
MAIN_DIR <- "C:/Users/alici/OneDrive/Documents/University/Year 2/internship!!!/Alicia/DATA"

# Define a folder for saving plots
PLOT_DIR <- paste0(MAIN_DIR, "/KM_CES1_plots")
if (!dir.exists(PLOT_DIR)) {
  dir.create(PLOT_DIR)
}

# Load COAD clinical dataset with CES1 expression
data <- read.csv(paste(MAIN_DIR, '/TOIL_RSEM_TPM_COAD_CES1_clinical.csv', sep = ""), 
                 sep = ",", 
                 row.names = NULL, 
                 header = TRUE)

dim(data)
head(data)

DF <- data
THRESHOLD <- 0.5 

# Create LOW/HIGH column
DF$CES1_level <- NA
DF$CES1_level[ DF$CES1_tpm <= quantile(DF$CES1_tpm, THRESHOLD, na.rm=TRUE) ] <- "LOW"
DF$CES1_level[ DF$CES1_tpm >  quantile(DF$CES1_tpm, THRESHOLD, na.rm=TRUE) ] <- "HIGH"

table(DF$CES1_level)

# === OS CURVE ===
sc <- "OS"
DF.os <- subset(DF, !is.na(OS) & !is.na(OS.time))
DF.os$CES1_level <- as.factor(DF.os[["CES1_level"]])

n_high <- sum(DF.os$CES1_level == "HIGH", na.rm = TRUE)
n_low  <- sum(DF.os$CES1_level == "LOW", na.rm = TRUE)

fit <- survfit(Surv(OS.time, OS) ~ CES1_level, data = DF.os)
p <- ggsurvplot(fit, xlab = "Days", ylab = "Overall Survival (OS)",
                legend.labs = c(paste0("HIGH (n = ", n_high, ")"),
                                paste0("LOW (n = ", n_low, ")")),
                legend.title = "", pval = TRUE, pval.coord = c(10, 0.1),
                pval.size = 6, size = 1.5, censor.size = 5)
print(p)

tiff(filename = paste(PLOT_DIR, "/KM_", sc, "_", THRESHOLD, "_v1.tiff", sep=""),
     width = 530, height = 500, units = "px", pointsize = 12, compression = "none", bg = "white")
print(p)
dev.off()

# === DSS CURVE ===
sc <- "DSS"
DF.dss <- subset(DF, !is.na(DSS) & !is.na(DSS.time))
DF.dss$CES1_level <- as.factor(DF.dss[["CES1_level"]])

n_high <- sum(DF.dss$CES1_level == "HIGH", na.rm = TRUE)
n_low  <- sum(DF.dss$CES1_level == "LOW", na.rm = TRUE)

fit <- survfit(Surv(DSS.time, DSS) ~ CES1_level, data = DF.dss)
p <- ggsurvplot(fit, xlab = "Days", ylab = "Disease-Specific Survival (DSS)",
                legend.labs = c(paste0("HIGH (n = ", n_high, ")"),
                                paste0("LOW (n = ", n_low, ")")),
                legend.title = "", pval = TRUE, pval.coord = c(10, 0.1),
                pval.size = 6, size = 1.5, censor.size = 5)
print(p)

tiff(filename = paste(PLOT_DIR, "/KM_", sc, "_", THRESHOLD, "_v1.tiff", sep=""),
     width = 530, height = 500, units = "px", pointsize = 12, compression = "none", bg = "white")
print(p)
dev.off()

# === PFI CURVE ===
sc <- "PFI"
DF.pfi <- subset(DF, !is.na(PFI) & !is.na(PFI.time))
DF.pfi$CES1_level <- as.factor(DF.pfi[["CES1_level"]])

n_high <- sum(DF.pfi$CES1_level == "HIGH", na.rm = TRUE)
n_low  <- sum(DF.pfi$CES1_level == "LOW", na.rm = TRUE)

fit <- survfit(Surv(PFI.time, PFI) ~ CES1_level, data = DF.pfi)
p <- ggsurvplot(fit, xlab = "Days", ylab = "Progression-Free Interval (PFI)",
                legend.labs = c(paste0("HIGH (n = ", n_high, ")"),
                                paste0("LOW (n = ", n_low, ")")),
                legend.title = "", pval = TRUE, pval.coord = c(10, 0.1),
                pval.size = 6, size = 1.5, censor.size = 5)
print(p)

tiff(filename = paste(PLOT_DIR, "/KM_", sc, "_", THRESHOLD, "_v1.tiff", sep=""),
     width = 530, height = 500, units = "px", pointsize = 12, compression = "none", bg = "white")
print(p)
dev.off()

# === SUBSET by Stage I/II/III ===
DF <- subset(data, !(ajcc_pathologic_tumor_stage %in% c("IVa", "IVb", "Stage IV", "Stage IVA", "Stage IVB", "Stage IVC", "", "[Discrepancy]", "[Unknown]", "Stage X")))
suffix <- "_pathStage123"
THRESHOLD <- 0.5

DF$CES1_level <- NA
DF$CES1_level[ DF$CES1_tpm <= quantile(DF$CES1_tpm, THRESHOLD, na.rm=TRUE) ] <- "LOW"
DF$CES1_level[ DF$CES1_tpm >  quantile(DF$CES1_tpm, THRESHOLD, na.rm=TRUE) ] <- "HIGH"

table(DF$CES1_level)

# OS
sc <- "OS"
DF.os <- subset(DF, !is.na(OS) & !is.na(OS.time))
DF.os$CES1_level <- as.factor(DF.os[["CES1_level"]])
n_high <- sum(DF.os$CES1_level == "HIGH", na.rm = TRUE)
n_low  <- sum(DF.os$CES1_level == "LOW", na.rm = TRUE)

fit <- survfit(Surv(OS.time, OS) ~ CES1_level, data = DF.os)
p <- ggsurvplot(fit, xlab = "Days", ylab = "Overall Survival (OS)",
                legend.labs = c(paste0("HIGH (n = ", n_high, ")"),
                                paste0("LOW (n = ", n_low, ")")),
                legend.title = "", pval = TRUE, pval.coord = c(10, 0.1),
                pval.size = 6, size = 1.5, censor.size = 5)
print(p)

tiff(filename = paste(PLOT_DIR, "/KM_", sc, "_", THRESHOLD, "_", suffix, "_v1.tiff", sep=""),
     width = 530, height = 500, units = "px", pointsize = 12, compression = "none", bg = "white")
print(p)
dev.off()

# DSS
sc <- "DSS"
DF.dss <- subset(DF, !is.na(DSS) & !is.na(DSS.time))
DF.dss$CES1_level <- as.factor(DF.dss[["CES1_level"]])
n_high <- sum(DF.dss$CES1_level == "HIGH", na.rm = TRUE)
n_low  <- sum(DF.dss$CES1_level == "LOW", na.rm = TRUE)

fit <- survfit(Surv(DSS.time, DSS) ~ CES1_level, data = DF.dss)
p <- ggsurvplot(fit, xlab = "Days", ylab = "Disease-Specific Survival (DSS)",
                legend.labs = c(paste0("HIGH (n = ", n_high, ")"),
                                paste0("LOW (n = ", n_low, ")")),
                legend.title = "", pval = TRUE, pval.coord = c(10, 0.1),
                pval.size = 6, size = 1.5, censor.size = 5)
print(p)

tiff(filename = paste(PLOT_DIR, "/KM_", sc, "_", THRESHOLD, "_", suffix, "_v1.tiff", sep=""),
     width = 530, height = 500, units = "px", pointsize = 12, compression = "none", bg = "white")
print(p)
dev.off()

# PFI
sc <- "PFI"
DF.pfi <- subset(DF, !is.na(PFI) & !is.na(PFI.time))
DF.pfi$CES1_level <- as.factor(DF.pfi[["CES1_level"]])
n_high <- sum(DF.pfi$CES1_level == "HIGH", na.rm = TRUE)
n_low  <- sum(DF.pfi$CES1_level == "LOW", na.rm = TRUE)

fit <- survfit(Surv(PFI.time, PFI) ~ CES1_level, data = DF.pfi)
p <- ggsurvplot(fit, xlab = "Days", ylab = "Progression-Free Interval (PFI)",
                legend.labs = c(paste0("HIGH (n = ", n_high, ")"),
                                paste0("LOW (n = ", n_low, ")")),
                legend.title = "", pval = TRUE, pval.coord = c(10, 0.1),
                pval.size = 6, size = 1.5, censor.size = 5)
print(p)

tiff(filename = paste(PLOT_DIR, "/KM_", sc, "_", THRESHOLD, "_", suffix, "_v1.tiff", sep=""),
     width = 530, height = 500, units = "px", pointsize = 12, compression = "none", bg = "white")
print(p)
dev.off()
