```R
# Import the data
load("Data/3T3.R")

## Time course of target genes of PPARg
par(mfcol=c(1,1))
Prediction <- Target_Genes$PPARG_1.motif
Mean <- log2(colMeans(2^GeneRPKM[ GeneRPKM$MinimalFDR <= 0.01 & GeneRPKM$Factor%in% Prediction[ Prediction$Target == 1, "Factor"],c(2:9)]))
SD <- apply(GeneRPKM[ GeneRPKM$MinimalFDR <= 0.01 & GeneRPKM$Factor%in% Prediction[ Prediction$Target == 1, "Factor"],c(2:9)],2,FUN="sd")
SD <- SD / sqrt(nrow(GeneRPKM[ GeneRPKM$MinimalFDR <= 0.01 & GeneRPKM$Factor%in% Prediction[ Prediction$Target == 1, "Factor"],]))
B <- barplot(Mean, las=1, main = "PPARG", ylim=c(0,8))
arrows(B, Mean - SD, B, Mean + SD, angle = 90, length = 0.05, code=3)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)