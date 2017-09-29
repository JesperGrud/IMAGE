```R
# Import the data
Pvalues <- read.delim("Data/Motif_Scoring_Pvalue.txt", header=T)

# Plot the ROC-like curve for ELF1
plot(Pvalues$ELF1_FPR, Pvalues$ELF1_TPR, type="o", pch=16, main="ELF1", ylab="True positive rate", xlab="False positive rate")
abline(0,1)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 2](../Links/Figure2.md)