```R
# Load the necessary libraries
library(edgeR)

# Import the IMAGE data
load("Data/hMSC.R")

# Import the count data
Counts <- read.delim("Data/Candidates.count.txt")
rownames(Counts) <- Counts$RefSeq

# Setup to analyze using edgeR
dd <- data.frame(a = gl(2,2))
Design <- model.matrix(~ a, dd)

# Analyze the count data using edgeR
DGE <- DGEList(Counts[,c(45,46,49,50)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
Control <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
Control$RefSeq <- rownames(Control)

DGE <- DGEList(Counts[,c(49,50,13,14)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
HSF <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
HSF$RefSeq <- rownames(HSF)

DGE <- DGEList(Counts[,c(49,50,19,20)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
MAZ <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
MAZ$RefSeq <- rownames(MAZ)

DGE <- DGEList(Counts[,c(49,50,25,26)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
MYBL1 <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
MYBL1$RefSeq <- rownames(MYBL1)

DGE <- DGEList(Counts[,c(49,50,37,38)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
SATB1 <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
SATB1$RefSeq <- rownames(SATB1)

DGE <- DGEList(Counts[,c(49,50,43,44)], group = c("TSHZ1","TSHZ1","Control","Control"))
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
TSHZ1 <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
TSHZ1$RefSeq <- rownames(TSHZ1)

## Setup to capture results
MeanDF <- data.frame(matrix(ncol=2, nrow=5))
SDDF <- data.frame(matrix(ncol=2, nrow=5))

## HSF1
# Setup up to run it
Targets <- Target_Genes$HSF1_1.motif
NumberTargets <- nrow(Targets[ Targets$Target == 1,])
NumberSignificant <- 1000

Genes <- HSF
Genes <- Genes[ order(Genes$FDR),]
Genes$Top <- 0
Genes[ c(1:1000),"Top"] <- 1

# Find the result using predictions
Prediction <- nrow(Targets[ Targets$Target == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])

# Randomize the targets, find overlap with significant genes
RandomTargets <- data.frame(matrix(ncol=1, nrow=1000))
for (m in 1:1000) {
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  RandomTargets[m,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])
}

# Randomize both groups
Randomized <- data.frame(matrix(ncol=1, nrow=1000))
for (q in 1:1000) {
  Genes$Random <- 0
  Genes[ sample(1:nrow(Genes),NumberSignificant),"Random"] <- 1
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  Randomized[q,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Random == 1,"RefSeq"],"Symbol"],])
}

# Save all results to a matrix
MeanDF[1,] <- c(mean((Prediction/1000)/(Randomized[,1]/1000)), mean((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))
SDDF[1,] <- c(sd((Prediction/1000)/(Randomized[,1]/1000)), sd((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))

## MYBL1
# Setup up to run it
Targets <- Target_Genes$MYBL1_1.motif
NumberTargets <- nrow(Targets[ Targets$Target == 1,])

Genes <- MYBL1
Genes <- Genes[ order(Genes$FDR),]
Genes$Top <- 0
Genes[ c(1:1000),"Top"] <- 1

# Find the result using predictions
Prediction <- nrow(Targets[ Targets$Target == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])

# Randomize the targets, find overlap with significant genes
RandomTargets <- data.frame(matrix(ncol=1, nrow=1000))
for (m in 1:1000) {
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  RandomTargets[m,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])
}

# Randomize both groups
Randomized <- data.frame(matrix(ncol=1, nrow=1000))
for (q in 1:1000) {
  Genes$Random <- 0
  Genes[ sample(1:nrow(Genes),NumberSignificant),"Random"] <- 1
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  Randomized[q,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Random == 1,"RefSeq"],"Symbol"],])
}

# Save all results to a matrix
MeanDF[2,] <- c(mean((Prediction/1000)/(Randomized[,1]/1000)), mean((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))
SDDF[2,] <- c(sd((Prediction/1000)/(Randomized[,1]/1000)), sd((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))

## MAZ
# Setup up to run it
Targets <- Target_Genes$MAZ_1.motif
NumberTargets <- nrow(Targets[ Targets$Target == 1,])

Genes <- MAZ
Genes <- Genes[ order(Genes$FDR),]
Genes$Top <- 0
Genes[ c(1:1000),"Top"] <- 1

# Find the result using predictions
Prediction <- nrow(Targets[ Targets$Target == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])

# Randomize the targets, find overlap with significant genes
RandomTargets <- data.frame(matrix(ncol=1, nrow=1000))
for (m in 1:1000) {
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  RandomTargets[m,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])
}

# Randomize both groups
Randomized <- data.frame(matrix(ncol=1, nrow=1000))
for (q in 1:1000) {
  Genes$Random <- 0
  Genes[ sample(1:nrow(Genes),NumberSignificant),"Random"] <- 1
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  Randomized[q,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Random == 1,"RefSeq"],"Symbol"],])
}

# Save all results to a matrix
MeanDF[3,] <- c(mean((Prediction/1000)/(Randomized[,1]/1000)), mean((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))
SDDF[3,] <- c(sd((Prediction/1000)/(Randomized[,1]/1000)), sd((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))

## SATB1
# Setup up to run it
Targets <- Target_Genes$ONECUT1_1.motif
NumberTargets <- nrow(Targets[ Targets$Target == 1,])

Genes <- SATB1
Genes <- Genes[ order(Genes$FDR),]
Genes$Top <- 0
Genes[ c(1:1000),"Top"] <- 1

# Find the result using predictions
Prediction <- nrow(Targets[ Targets$Target == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])

# Randomize the targets, find overlap with significant genes
RandomTargets <- data.frame(matrix(ncol=1, nrow=1000))
for (m in 1:1000) {
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  RandomTargets[m,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])
}

# Randomize both groups
Randomized <- data.frame(matrix(ncol=1, nrow=1000))
for (q in 1:1000) {
  Genes$Random <- 0
  Genes[ sample(1:nrow(Genes),NumberSignificant),"Random"] <- 1
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  Randomized[q,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Random == 1,"RefSeq"],"Symbol"],])
}

# Save all results to a matrix
MeanDF[4,] <- c(mean((Prediction/1000)/(Randomized[,1]/1000)), mean((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))
SDDF[4,] <- c(sd((Prediction/1000)/(Randomized[,1]/1000)), sd((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))

## TSHZ1
# Setup up to run it
Targets <- Target_Genes$ZEB1_1.motif
NumberTargets <- nrow(Targets[ Targets$Target == 1,])

Genes <- TSHZ1
Genes <- Genes[ order(Genes$FDR),]
Genes$Top <- 0
Genes[ c(1:1000),"Top"] <- 1

# Find the result using predictions
Prediction <- nrow(Targets[ Targets$Target == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])

# Randomize the targets, find overlap with significant genes
RandomTargets <- data.frame(matrix(ncol=1, nrow=1000))
for (m in 1:1000) {
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  RandomTargets[m,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Top == 1,"RefSeq"],"Symbol"],])
}

# Randomize both groups
Randomized <- data.frame(matrix(ncol=1, nrow=1000))
for (q in 1:1000) {
  Genes$Random <- 0
  Genes[ sample(1:nrow(Genes),NumberSignificant),"Random"] <- 1
  Targets$Random <- 0
  Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
  Randomized[q,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Counts[ Counts$RefSeq %in% Genes[ Genes$Random == 1,"RefSeq"],"Symbol"],])
}

# Save all results to a matrix
MeanDF[5,] <- c(mean((Prediction/1000)/(Randomized[,1]/1000)), mean((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))
SDDF[5,] <- c(sd((Prediction/1000)/(Randomized[,1]/1000)), sd((mean(RandomTargets[,1]/1000))/(Randomized[,1]/1000)))

## Plot the results
# Setup
MeanDF <- t(as.matrix(MeanDF))
SDDF <- t(as.matrix(SDDF))

# Plot it
par(mfcol=c(1,1))
B <- barplot(log2(MeanDF), beside=T, names=c("HSF1","MYBL1","MAZ","SATB1","TSHZ1"), las=2, ylab="log2 Enrichment", ylim=c(0,2.5), col=c("darkblue","lightgrey"))
arrows(B, log2(MeanDF+SDDF), B, log2(MeanDF-SDDF), code=3, angle=90, length=0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure S2](../Links/FigureS2.md)