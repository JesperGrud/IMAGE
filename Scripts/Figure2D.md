```R
# Read the data
Pvalues <- read.delim("Data/Motif_Scoring_Pvalue.txt", header=T)

# Generate sensitivity times specificity matrix
SensSpec <- data.frame(matrix(ncol=15,nrow=20))
SensSpec[,1] <- Pvalues[,1]
Counter <- 2

for (i in seq(2,ncol(Pvalues),by=2)) {
SensSpec[,Counter] <- Pvalues[,i]*(1-Pvalues[,(i+1)])
Counter <- Counter + 1
}

# Plot the results
plot(-log10(SensSpec[,1]), rowMeans(SensSpec[,c(2:15)]), type="l", ylim=c(0,1), xlab="-log10(Pvalue)", ylab = "Sensitivity * Specificity", col="green", lwd=2)
arrows(-log10(SensSpec[,1]), rowMeans(SensSpec[,c(2:15)])-apply(SensSpec[,c(2:15)],1,FUN="sd"),-log10(SensSpec[,1]), rowMeans(SensSpec[,c(2:15)])+apply(SensSpec[,c(2:15)],1,FUN="sd"), code=3, angle=90, length = 0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 2](../Links/Figure2.md)