```R
# Source the necessary libraries
library(beeswarm)

# Read the data
Pvalues <- read.delim("Data/Motif_Scoring_Pvalue.txt", header=T)
LogRT <- read.delim("Data/Motif_Scoring_LogRatio.txt", header=T)

# Generate sensitivity times specificity matrix for p-value-based scoring
SensSpecPvalues <- data.frame(matrix(ncol=15,nrow=20))
SensSpecPvalues[,1] <- Pvalues[,1]
Counter <- 2

for (i in seq(2,ncol(Pvalues),by=2)) {
SensSpecPvalues[,Counter] <- Pvalues[,i]*(1-Pvalues[,(i+1)])
Counter <- Counter + 1
}

# Generate sensitivity times specificity matrix for log-likelihood scoring
SensSpecLL <- data.frame(matrix(ncol=15,nrow=35))
SensSpecLL[,1] <- LogRT[,1]
Counter <- 2

for (i in seq(2,ncol(LogRT),by=2)) {
  SensSpecLL[,Counter] <- LogRT[,i]*(1-LogRT[,(i+1)])
  Counter <- Counter + 1
}

# For each motif get the relative performance at the common max
RelativePerformance <- data.frame(matrix(ncol=2, nrow=14))
RelativePerformance[,1] <- t(SensSpecPvalues[which.max(rowMeans(SensSpecPvalues[,c(2:15)])),c(2:15)]/apply(SensSpecPvalues[,c(2:15)],2,FUN="max"))
RelativePerformance[,2] <- t(SensSpecLL[which.max(rowMeans(SensSpecLL[,c(2:15)])),c(2:15)]/apply(SensSpecLL[,c(2:15)],2,FUN="max"))
colnames(RelativePerformance) <- c("P-value","Log-likelihood")

# Plot the results
beeswarm(RelativePerformance, ylab="Relative performance", pch=16, cex=1, col="black", las=1)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 2](../Links/Figure2.md)