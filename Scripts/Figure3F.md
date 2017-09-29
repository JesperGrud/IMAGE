```R
# Import the data
load("Data/3T3.R")

## Enhancer time course for PPARg predicted sites
Enhancers$D0 <- rowMeans(Enhancers[,c(5,6)])
Enhancers$H4 <- rowMeans(Enhancers[,c(7,8)])
Enhancers$D2 <- rowMeans(Enhancers[,c(9,10)])
Enhancers$D7 <- rowMeans(Enhancers[,c(11,12)])

par(mfcol=c(1,2))
All <- Target_Sites$PPARG_1.motif
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
boxplot(Enhancers[ Enhancers$ID %in% Model$PeakID, c(13,14,15,16)], outline=F, las=1, names=c("D0","4h","D2","D7"))
All <- Target_Sites$NR3C1_1.motif
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
boxplot(Enhancers[ Enhancers$ID %in% Model$PeakID, c(13,14,15,16)], outline=F, las=1, names=c("D0","4h","D2","D7"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 3](../Links/Figure3.md)