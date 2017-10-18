```R
# Import the data
load("Data/3T3.R")

# Import some data
MotifToFactor <- read.delim("Data/Genename_Motif.txt", header=FALSE)
TFs <- read.table("Data/TFs.txt", quote="\"", comment.char="")

# All expressed
Expressed <- data.frame(Factor = TFs[ TFs[,1] %in% GeneRPKM[ GeneRPKM$Max >= 0,"Factor"],])
Expressed <- data.frame(Factor = Expressed[ duplicated(Expressed$Factor)==F,])

# Regulated during differentiation
ExpressedRegulated <- data.frame(Factor = TFs[ TFs[,1] %in% GeneRPKM[ GeneRPKM$Max >= 0 & ((GeneRPKM$`FDR_Cond3-vs-Cond2` <= 0.05) | (GeneRPKM$`FDR_Cond4-vs-Cond2` <= 0.05) | (GeneRPKM$`FDR_Cond4-vs-Cond3` <= 0.05) | (GeneRPKM$`FDR_Cond4-vs-Cond1` <= 0.05)  | (GeneRPKM$`FDR_Cond3-vs-Cond1` <= 0.05)  | (GeneRPKM$`FDR_Cond2-vs-Cond1` <= 0.05)),"Factor"],])

## Make groups of transcription factors based on motif enrichments
MotifsAll <- read.delim("Data/All_knownResults.txt")
MotifsAll <- MotifsAll[ MotifsAll[,3] <= 0.05, ]
MotifsAll <- merge(MotifsAll, MotifToFactor, by.x="Motif.Name", by.y="V2")
MotifsAll <- MotifsAll[ duplicated(MotifsAll$V1)==F,]

MotifsInd <- read.delim("Data/Induced_knownResults.txt")
MotifsInd <- MotifsInd[ MotifsInd[,3] <= 0.05, c(1:5)]
MotifsRep <- read.delim("Data/Repressed_knownResults.txt")
MotifsRep <- MotifsRep[ MotifsRep[,3] <= 0.05, c(1:5)]
MotifsRegulated <- rbind(MotifsInd,MotifsRep)
MotifsRegulated <- merge(MotifsRegulated, MotifToFactor, by.x="Motif.Name", by.y="V2")
MotifsRegulated <- MotifsRegulated[ duplicated(MotifsRegulated$V1)==F,]

## Make groups based on combinations of motifs and expression
Combined <- MotifsRegulated[ MotifsRegulated$V1 %in% ExpressedRegulated$Factor,]

## Process the expression data
Data <- read.delim("Data/GDS3961_full.soft")
Samples <- read.delim("Data/AdiposeTissueSamples.txt")
Data <- Data[ ,colnames(Data) %in% Samples[,1] | colnames(Data) == "IDENTIFIER" | colnames(Data) == "ID_REF"]

lmEC <- lm(Samples[,"BMI"] ~ Samples[,"Age"])
lmEC.resid <- residuals(lmEC)

GOI <- data.frame(Symbol = TFs$V1)
GOI <- data.frame(Symbol = GOI[ GOI$Symbol %in% Data$IDENTIFIER,])
Enrichment <- data.frame()

for (i in 1:nrow(GOI)) {
  NCol <- 1
  Gene <- GOI[i,1]
  Correlate <- as.data.frame(t(Data[ Data$IDENTIFIER == as.character(Gene), ]))
  Correlate$Sample <- rownames(Correlate)
  Correlate <- Correlate[c(2:nrow(Correlate)),]
  for (i in 1:(ncol(Correlate)-1)) { Correlate[,i] <- as.numeric(as.character(Correlate[,i]))}
  if (ncol(Correlate) > 2) {
    NCol <- ncol(Correlate)
    Correlate$Mean <- apply(Correlate[,c(1:(ncol(Correlate)-1))],1,FUN="mean")
  }
  Correlate <- merge(Correlate, Samples[,c(1,3,4,2)], by="Sample")
  
  Tmp <- data.frame()
  for (i in 1:NCol) {
  Tmp[i,1] <- summary(lm(Correlate[,(i+1)] ~ lmEC.resid))$r.square
  Tmp[i,2] <- summary(lm(Correlate[,(i+1)] ~ lmEC.resid))$coefficients[2,4]
  }
  Tmp$Row <- which.min(Tmp[,2])
  Tmp <- Tmp[ which.min(Tmp$V2),]
  Tmp$Gene <- Gene
  Enrichment <- rbind(Enrichment, Tmp)
  rm(Tmp)
}

# Calculate proportions
Matrix <- data.frame(matrix(ncol=3, nrow=8))
Matrix[1,2] <- nrow(Enrichment[ Enrichment$Gene %in% Result[ Result$CausalTF == 1, "Factor"],])
Matrix[1,1] <- nrow(Enrichment[ Enrichment$Gene %in% Result[ Result$CausalTF == 1, "Factor"] & Enrichment[,2] <= 0.05,])
Matrix[2,2] <- nrow(Enrichment[ Enrichment$Gene %in% Combined$V1,])
Matrix[2,1] <- nrow(Enrichment[ Enrichment$Gene %in% Combined$V1 & Enrichment[,2] <= 0.05,])
Matrix[3,2] <- nrow(Enrichment[ Enrichment$Gene %in% MotifsRegulated$V1,])
Matrix[3,1] <- nrow(Enrichment[ Enrichment$Gene %in% MotifsRegulated$V1 & Enrichment[,2] <= 0.05,])
Matrix[4,2] <- nrow(Enrichment[ Enrichment$Gene %in% MotifsAll$V1,])
Matrix[4,1] <- nrow(Enrichment[ Enrichment$Gene %in% MotifsAll$V1 & Enrichment[,2] <= 0.05,])
Matrix[5,2] <- nrow(Enrichment[ Enrichment$Gene %in% ExpressedRegulated$Factor,])
Matrix[5,1] <- nrow(Enrichment[ Enrichment$Gene %in% ExpressedRegulated$Factor & Enrichment[,2] <= 0.05,])
Matrix[6,2] <- nrow(Enrichment[ Enrichment$Gene %in% Expressed$Factor,])
Matrix[6,1] <- nrow(Enrichment[ Enrichment$Gene %in% Expressed$Factor & Enrichment[,2] <= 0.05,])
Matrix[8,2] <- nrow(Enrichment[ Enrichment$Gene %in% TFs$V1,])
Matrix[8,1] <- nrow(Enrichment[ Enrichment$Gene %in% TFs$V1 & Enrichment[,2] <= 0.05,])

# Perform permutations on random selection
RandomResult <- data.frame(matrix(ncol=2, nrow=1000))
for (i in 1:1000) {
  Random <- data.frame(Factor = TFs[ sample(1:nrow(TFs),500),1])
  RandomResult[i,1] <- nrow(Enrichment[ Enrichment$Gene %in% Random$Factor & Enrichment[,2] <= 0.05,])
  RandomResult[i,2] <- nrow(Enrichment[ Enrichment$Gene %in% Random$Factor,])
}

Matrix[7,1] <- mean(RandomResult[,1])
Matrix[7,2] <- mean(RandomResult[,2])

# Calculate enrichments
Matrix[,3] <- Matrix[,1] / Matrix[,2]
Matrix[,3] <- Matrix[,3] / Matrix[8,3]
Matrix[,3] <- log2(Matrix[,3])

# Make the plot
par(mfcol=c(1,1))
barplot(Matrix[c(7:1),3], las=2, ylim=c(0,1.2), ylab="Log2 Enrichment over all TFs", names=c("Random","All","Dynamic","All","Dynamic","dMotif+cExprs","IMAGE"), col=c("grey","blue","blue","red","red","purple","green"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)