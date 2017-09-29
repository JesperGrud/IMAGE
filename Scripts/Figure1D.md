```R
# Import the data
Correlation <- read.delim("Data/Prediction_Validation.txt", header=T)
Correlation <- Correlation[ rank(-Correlation$Cor),]

# Plot the results
plot(Correlation$Cor, type="l", las=1, ylim=c(0,1), col="red", lwd=2, xlab="Number", ylab="Pearson's correlation coefficient")
abline(h=0.8, lty=2)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 1](../Links/Figure1.md)