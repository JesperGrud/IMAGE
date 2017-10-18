```R
# Import the data
load("Data/3T3.R")
MED1 <- Target_Sites
rm(list=setdiff(ls(),c("MED1")))
load("Data/DHS_Filtered.RData")
DHS <- Target_Sites
rm(Target_Sites)

# Get motifs checked in both data sets
MED1Names <- data.frame(Motif = names(MED1))
DHSNames <- data.frame(Motif = names(DHS))
Motifs <- data.frame(Motif = as.character(MED1Names[ MED1Names$Motif %in% DHSNames$Motif,]))
rm(MED1Names)
rm(DHSNames)

# Loop through all motifs and calculate overlap for predicted sites
Result <- as.data.frame(matrix(ncol=4, nrow=nrow(Motifs)))
colnames(Result) <- c("Motif","Shared","DHS_unique","MED1_unique")
for (i in 1:nrow(Motifs)) {
Motif <- as.character(Motifs[i,1])
MED1Sites <- as.data.frame(MED1[names(MED1) == Motif])
MED1Sites <- MED1Sites[ MED1Sites[,5] > 0 & MED1Sites[,6] > 0,]
DHSSites <- as.data.frame(DHS[names(DHS) == Motif])
DHSSites <- DHSSites[ DHSSites[,5] > 0 & DHSSites[,6] > 0,]
Result[i,1] <- Motif
Result[i,2] <- nrow(DHSSites[ DHSSites[,4] %in% MED1Sites[,4], ])
Result[i,3] <- nrow(DHSSites[ !(DHSSites[,4] %in% MED1Sites[,4]), ])
Result[i,4] <- nrow(MED1Sites[ !(MED1Sites[,4] %in% DHSSites[,4]), ])
}

Result$Sum <- apply(Result[,c(2:4)],1,FUN="sum")

# Predict randomly
Result$MED1 <- Result[,2] + Result[,4]
Result$DHS <- Result[,2] + Result[,5]
Sites <- DHS[[1]]

for (i in 1:nrow(Result)) {
Tmp <- Sites[ sample(1:nrow(Sites), Result[i,"MED1"]),]
Tmp2 <- Sites[ sample(1:nrow(Sites), Result[i,"DHS"]),]
Result[i,8] <- nrow(Tmp2[ Tmp2[,4] %in% Tmp[,4], ])
Result[i,9] <- nrow(Tmp2[ !(Tmp2[,4] %in% Tmp[,4]), ])
Result[i,10] <- nrow(Tmp[ !(Tmp[,4] %in% Tmp2[,4]), ])
}
Result$Sum2 <- apply(Result[,c(8:10)],1,FUN="sum")

# Plot the results
boxplot(Result[,8]/Result[,11], Result[,2]/Result[,5], outline=F, col="grey", las=1, ylab="Fraction", names=c("Random","Predicted"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)