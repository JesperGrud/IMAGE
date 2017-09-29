```R
# Read the data
LogRT <- read.delim("Data/Motif_Scoring_LogRatio.txt", header=T)

# Generate sensitivity times specificity matrix
SensSpec <- data.frame(matrix(ncol=15,nrow=35))
SensSpec[,1] <- LogRT[,1]
Counter <- 2

for (i in seq(2,ncol(LogRT),by=2)) {
  SensSpec[,Counter] <- LogRT[,i]*(1-LogRT[,(i+1)])
  Counter <- Counter + 1
}

# Plot the results
plot(SensSpec[,1], rowMeans(SensSpec[,c(2:15)]), type="l", ylim=c(-0.1,1), xlab="log-likelihood", ylab = "Sensitivity * Specificity", col="green", lwd=2)
arrows(SensSpec[,1], rowMeans(SensSpec[,c(2:15)])-apply(SensSpec[,c(2:15)],1,FUN="sd"),SensSpec[,1], rowMeans(SensSpec[,c(2:15)])+apply(SensSpec[,c(2:15)],1,FUN="sd"), code=3, angle=90, length = 0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)