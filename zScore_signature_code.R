###survival code
library(survival)
library(survminer)

# 4a: Signature genes
SIGNATUselected.tmp <- DF[, SIGNATURE]RE <- c("CES1", "CPT1A")


# 4a: Extract and coerce safely
selected.tmp <- subset(DF, select = intersect(SIGNATURE, names(DF)))
selected.tmp[] <- lapply(selected.tmp, function(x) as.numeric(as.character(x)))

# Drop rows with any NA in signature genes
valid_rows <- complete.cases(selected.tmp)
selected.tmp <- selected.tmp[valid_rows, ]
DF <- DF[valid_rows, ]

# 4a: Z-score normalize + compute composite
myDF.zscore <- scale(selected.tmp)
mean.zscore <- rowMeans(myDF.zscore)
DF$zsign_CES1_CPT1A <- mean.zscore

# 4b: Create thresholds
DF$zsign_CES1_CPT1A_25p <- ifelse(
  DF$zsign_CES1_CPT1A > quantile(DF$zsign_CES1_CPT1A, 0.25, na.rm = TRUE),
  "HIGH", "LOW"
)

# 4c: Survival analysis
gg <- "zsign_CES1_CPT1A_25p"
thres <- 0
tmp.sub <- subset(DF, !is.na(OS) & !is.na(OS.time) & OS.time >= thres)

tmp.sub$usedSubset <- as.factor(tmp.sub[[gg]])
tmp.sub <- subset(tmp.sub, !is.na(usedSubset))

fit <- survfit(Surv(OS.time, OS) ~ usedSubset, data = tmp.sub)

# Dynamic legend
groups_present <- levels(droplevels(tmp.sub$usedSubset))
legend_labels <- c()
if ("HIGH" %in% groups_present) legend_labels <- c(legend_labels, paste0("High (n = ", sum(tmp.sub$usedSubset == "HIGH"), ")"))
if ("LOW" %in% groups_present)  legend_labels <- c(legend_labels, paste0("Low (n = ", sum(tmp.sub$usedSubset == "LOW"), ")"))

p <- ggsurvplot(
  fit,
  palette = c("#ED0000FF", "black"),
  title = paste0("Survival by ", gg),
  xlab = "Months",
  ylab = "Overall Survival (OS)",
  break.time.by = 30,
  legend = c(0.6, 0.9),
  legend.title = "",
  legend.labs = legend_labels,
  pval = TRUE,
  pval.size = 5,
  fun = function(y) y * 100
)

print(p)


CES1 -> ENSG00000198848.12
CPT1A -> ENSG00000110090.12
MGLL -> ENSG00000074416.13
ACSL3 -> ENSG00000123983.13


SIGNATURE <- c("CES1","CPT1A")
selected.tmp <- subset(DF, select=c(intersect(SIGNATURE, names(DF))) ) #  
myDF.zscore <- scale( selected.tmp)
mean.zscore <- as.data.frame(apply(myDF.zscore, 1, mean  ) )
names(mean.zscore) <- c( "zsign_CES1_CPT1A")

DF <- merge(subset(DF ),mean.zscore,by="row.names")
row.names(DF) <- DF$Sample 
DF <- subset(DF, select=-c(Row.names))

DF$zsign_CES1_CPT1A_25p<- "LOW"
DF$zsign_CES1_CPT1A_25p[DF$zsign_CES1_CPT1A  > quantile( DF$zsign_CES1_CPT1A,0.25) ]<- "HIGH"
DF$zsign_CES1_CPT1A_50p<- "LOW"
DF$zsign_CES1_CPT1A_50p[DF$zsign_CES1_CPT1A  > quantile( DF$zsign_CES1_CPT1A,0.50) ]<- "HIGH"
DF$zsign_CES1_CPT1A_75p<- "LOW"
DF$zsign_CES1_CPT1A_75p[DF$zsign_CES1_CPT1A  > quantile( DF$zsign_CES1_CPT1A,0.75) ]<- "HIGH"
DF$zsign_CES1_CPT1A_66p<- "LOW"
DF$zsign_CES1_CPT1A_66p[DF$zsign_CES1_CPT1A  > quantile( DF$zsign_CES1_CPT1A,0.66) ]<- "HIGH"
DF$zsign_CES1_CPT1A_90p<- "LOW"
DF$zsign_CES1_CPT1A_90p[DF$zsign_CES1_CPT1A  > quantile( DF$zsign_CES1_CPT1A,0.90) ]<- "HIGH"

DF$zsign_CES1_CPT1A_25p_75p<- NA
DF$zsign_CES1_CPT1A_25p_75p[DF$zsign_CES1_CPT1A  < quantile( DF$zsign_CES1_CPT1A,0.25) ]<- "LOW"
DF$zsign_CES1_CPT1A_25p_75p[DF$zsign_CES1_CPT1A  > quantile( DF$zsign_CES1_CPT1A,0.75) ]<- "HIGH"  
DF$zsign_CES1_CPT1A_33p_66p<- NA
DF$zsign_CES1_CPT1A_33p_66p[DF$zsign_CES1_CPT1A  < quantile( DF$zsign_CES1_CPT1A,0.33) ]<- "LOW"
DF$zsign_CES1_CPT1A_33p_66p[DF$zsign_CES1_CPT1A  > quantile( DF$zsign_CES1_CPT1A,0.66) ]<- "HIGH" 

SIGNATURE <- c("CES1","MGLL","ACSL3")
selected.tmp <- subset(DF, select=c(intersect(SIGNATURE, names(DF))) ) #  
myDF.zscore <- scale( selected.tmp)
mean.zscore <- as.data.frame(apply(myDF.zscore, 1, mean  ) )
names(mean.zscore) <- c( "zsign_CES1_MGLL_ACSL3")

DF <- merge(subset(DF ),mean.zscore,by="row.names")
row.names(DF) <- DF$Sample 
DF <- subset(DF, select=-c(Row.names))

DF$zsign_CES1_MGLL_ACSL3_25p<- "LOW"
DF$zsign_CES1_MGLL_ACSL3_25p[DF$zsign_CES1_MGLL_ACSL3  > quantile( DF$zsign_CES1_MGLL_ACSL3,0.25) ]<- "HIGH"
DF$zsign_CES1_MGLL_ACSL3_50p<- "LOW"
DF$zsign_CES1_MGLL_ACSL3_50p[DF$zsign_CES1_MGLL_ACSL3  > quantile( DF$zsign_CES1_MGLL_ACSL3,0.50) ]<- "HIGH"
DF$zsign_CES1_MGLL_ACSL3_75p<- "LOW"
DF$zsign_CES1_MGLL_ACSL3_75p[DF$zsign_CES1_MGLL_ACSL3  > quantile( DF$zsign_CES1_MGLL_ACSL3,0.75) ]<- "HIGH"
DF$zsign_CES1_MGLL_ACSL3_66p<- "LOW"
DF$zsign_CES1_MGLL_ACSL3_66p[DF$zsign_CES1_MGLL_ACSL3  > quantile( DF$zsign_CES1_MGLL_ACSL3,0.66) ]<- "HIGH"
DF$zsign_CES1_MGLL_ACSL3_90p<- "LOW"
DF$zsign_CES1_MGLL_ACSL3_90p[DF$zsign_CES1_MGLL_ACSL3  > quantile( DF$zsign_CES1_MGLL_ACSL3,0.90) ]<- "HIGH"

DF$zsign_CES1_MGLL_ACSL3_25p_75p<- NA
DF$zsign_CES1_MGLL_ACSL3_25p_75p[DF$zsign_CES1_MGLL_ACSL3  < quantile( DF$zsign_CES1_MGLL_ACSL3,0.25) ]<- "LOW"
DF$zsign_CES1_MGLL_ACSL3_25p_75p[DF$zsign_CES1_MGLL_ACSL3  > quantile( DF$zsign_CES1_MGLL_ACSL3,0.75) ]<- "HIGH"  
DF$zsign_CES1_MGLL_ACSL3_33p_66p<- NA
DF$zsign_CES1_MGLL_ACSL3_33p_66p[DF$zsign_CES1_MGLL_ACSL3  < quantile( DF$zsign_CES1_MGLL_ACSL3,0.33) ]<- "LOW"
DF$zsign_CES1_MGLL_ACSL3_33p_66p[DF$zsign_CES1_MGLL_ACSL3  > quantile( DF$zsign_CES1_MGLL_ACSL3,0.66) ]<- "HIGH"  

SIGNATURE <- c("CES1","MGLL","CPT1A")
selected.tmp <- subset(DF, select=c(intersect(SIGNATURE, names(DF))) ) #  
myDF.zscore <- scale( selected.tmp)
mean.zscore <- as.data.frame(apply(myDF.zscore, 1, mean  ) )
names(mean.zscore) <- c( "zsign_CES1_MGLL_CPT1A")

DF <- merge(subset(DF ),mean.zscore,by="row.names")
row.names(DF) <- DF$Sample 
DF <- subset(DF, select=-c(Row.names))

DF$zsign_CES1_MGLL_CPT1A_25p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_25p[DF$zsign_CES1_MGLL_CPT1A  > quantile( DF$zsign_CES1_MGLL_CPT1A,0.25) ]<- "HIGH"
DF$zsign_CES1_MGLL_CPT1A_50p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_50p[DF$zsign_CES1_MGLL_CPT1A  > quantile( DF$zsign_CES1_MGLL_CPT1A,0.50) ]<- "HIGH"
DF$zsign_CES1_MGLL_CPT1A_75p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_75p[DF$zsign_CES1_MGLL_CPT1A  > quantile( DF$zsign_CES1_MGLL_CPT1A,0.75) ]<- "HIGH"
DF$zsign_CES1_MGLL_CPT1A_66p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_66p[DF$zsign_CES1_MGLL_CPT1A  > quantile( DF$zsign_CES1_MGLL_CPT1A,0.66) ]<- "HIGH"
DF$zsign_CES1_MGLL_CPT1A_90p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_90p[DF$zsign_CES1_MGLL_CPT1A  > quantile( DF$zsign_CES1_MGLL_CPT1A,0.90) ]<- "HIGH"

DF$zsign_CES1_MGLL_CPT1A_25p_75p<- NA
DF$zsign_CES1_MGLL_CPT1A_25p_75p[DF$zsign_CES1_MGLL_CPT1A  < quantile( DF$zsign_CES1_MGLL_CPT1A,0.25) ]<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_25p_75p[DF$zsign_CES1_MGLL_CPT1A  > quantile( DF$zsign_CES1_MGLL_CPT1A,0.75) ]<- "HIGH"  
DF$zsign_CES1_MGLL_CPT1A_33p_66p<- NA
DF$zsign_CES1_MGLL_CPT1A_33p_66p[DF$zsign_CES1_MGLL_CPT1A  < quantile( DF$zsign_CES1_MGLL_CPT1A,0.33) ]<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_33p_66p[DF$zsign_CES1_MGLL_CPT1A  > quantile( DF$zsign_CES1_MGLL_CPT1A,0.66) ]<- "HIGH"  

SIGNATURE <- c("CES1","MGLL","CPT1A", "ACSL3")
selected.tmp <- subset(DF, select=c(intersect(SIGNATURE, names(DF))) ) #  
myDF.zscore <- scale( selected.tmp)
mean.zscore <- as.data.frame(apply(myDF.zscore, 1, mean  ) )
names(mean.zscore) <- c( "zsign_CES1_MGLL_CPT1A_ACSL3")

DF <- merge(subset(DF ),mean.zscore,by="row.names")
row.names(DF) <- DF$Sample 
DF <- subset(DF, select=-c(Row.names))

DF$zsign_CES1_MGLL_CPT1A_ACSL3_25p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_25p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  > quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.25) ]<- "HIGH"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_50p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_50p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  > quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.50) ]<- "HIGH"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_75p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_75p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  > quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.75) ]<- "HIGH"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_66p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_66p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  > quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.66) ]<- "HIGH"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_90p<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_90p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  > quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.90) ]<- "HIGH"

DF$zsign_CES1_MGLL_CPT1A_ACSL3_25p_75p<- NA
DF$zsign_CES1_MGLL_CPT1A_ACSL3_25p_75p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  < quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.25) ]<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_25p_75p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  > quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.75) ]<- "HIGH"  
DF$zsign_CES1_MGLL_CPT1A_ACSL3_33p_66p<- NA
DF$zsign_CES1_MGLL_CPT1A_ACSL3_33p_66p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  < quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.33) ]<- "LOW"
DF$zsign_CES1_MGLL_CPT1A_ACSL3_33p_66p[DF$zsign_CES1_MGLL_CPT1A_ACSL3  > quantile( DF$zsign_CES1_MGLL_CPT1A_ACSL3,0.66) ]<- "HIGH"  

SIGNATURE <- c("CES1", "CPT1A", "ACSL3")
selected.tmp <- subset(DF, select=c(intersect(SIGNATURE, names(DF))) ) #  
myDF.zscore <- scale( selected.tmp)
mean.zscore <- as.data.frame(apply(myDF.zscore, 1, mean  ) )
names(mean.zscore) <- c( "zsign_CES1_CPT1A_ACSL3")

DF <- merge(subset(DF ),mean.zscore,by="row.names")
row.names(DF) <- DF$Sample 
DF <- subset(DF, select=-c(Row.names))

DF$zsign_CES1_CPT1A_ACSL3_25p<- "LOW"
DF$zsign_CES1_CPT1A_ACSL3_25p[DF$zsign_CES1_CPT1A_ACSL3  > quantile( DF$zsign_CES1_CPT1A_ACSL3,0.25) ]<- "HIGH"
DF$zsign_CES1_CPT1A_ACSL3_50p<- "LOW"
DF$zsign_CES1_CPT1A_ACSL3_50p[DF$zsign_CES1_CPT1A_ACSL3  > quantile( DF$zsign_CES1_CPT1A_ACSL3,0.50) ]<- "HIGH"
DF$zsign_CES1_CPT1A_ACSL3_75p<- "LOW"
DF$zsign_CES1_CPT1A_ACSL3_75p[DF$zsign_CES1_CPT1A_ACSL3  > quantile( DF$zsign_CES1_CPT1A_ACSL3,0.75) ]<- "HIGH"
DF$zsign_CES1_CPT1A_ACSL3_66p<- "LOW"
DF$zsign_CES1_CPT1A_ACSL3_66p[DF$zsign_CES1_CPT1A_ACSL3  > quantile( DF$zsign_CES1_CPT1A_ACSL3,0.66) ]<- "HIGH"
DF$zsign_CES1_CPT1A_ACSL3_90p<- "LOW"
DF$zsign_CES1_CPT1A_ACSL3_90p[DF$zsign_CES1_CPT1A_ACSL3  > quantile( DF$zsign_CES1_CPT1A_ACSL3,0.90) ]<- "HIGH"

DF$zsign_CES1_CPT1A_ACSL3_25p_75p<- NA
DF$zsign_CES1_CPT1A_ACSL3_25p_75p[DF$zsign_CES1_CPT1A_ACSL3  < quantile( DF$zsign_CES1_CPT1A_ACSL3,0.25) ]<- "LOW"
DF$zsign_CES1_CPT1A_ACSL3_25p_75p[DF$zsign_CES1_CPT1A_ACSL3  > quantile( DF$zsign_CES1_CPT1A_ACSL3,0.75) ]<- "HIGH"  
DF$zsign_CES1_CPT1A_ACSL3_33p_66p<- NA
DF$zsign_CES1_CPT1A_ACSL3_33p_66p[DF$zsign_CES1_CPT1A_ACSL3  < quantile( DF$zsign_CES1_CPT1A_ACSL3,0.33) ]<- "LOW"
DF$zsign_CES1_CPT1A_ACSL3_33p_66p[DF$zsign_CES1_CPT1A_ACSL3  > quantile( DF$zsign_CES1_CPT1A_ACSL3,0.66) ]<- "HIGH"  