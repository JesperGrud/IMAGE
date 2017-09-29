```R
# Import the data
load("Data/3T3.R")

# Import the list of GO terms and transcription factors
GO <- read.delim("Data/FatCell_GO.txt", header=TRUE, stringsAsFactors=FALSE)
TFs <- read.table("Data/TFs.txt", quote="\"", comment.char="")
MotifToFactor <- read.delim("Data/Genename_Motif.txt", header=FALSE)

## Make groups of transcription factors based on expression patterns
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

## Collaps hits
Hits <- Result[ Result$CausalTF == 1,]
Hits <- Hits[ duplicated(Hits$Factor)==F,]

# Calculate enrichment for different groups and store in a matrix
Enrichment <- data.frame(Fraction = matrix(ncol=2,nrow=7))
Enrichment[1,1] <- nrow(Hits[ Hits$Factor %in% GO$Factor,])
Enrichment[1,2] <- nrow(Hits)
Enrichment[2,1] <- nrow(Combined[ Combined$V1 %in% GO$Factor,])
Enrichment[2,2] <- nrow(Combined)
Enrichment[3,1] <- nrow(MotifsRegulated[ MotifsRegulated$V1 %in% GO$Factor,])
Enrichment[3,2] <- nrow(MotifsRegulated)
Enrichment[4,1] <- nrow(MotifsAll[ MotifsAll$V1 %in% GO$Factor,])
Enrichment[4,2] <- nrow(MotifsAll)
Enrichment[5,1] <- length(ExpressedRegulated[ ExpressedRegulated$Factor %in% GO$Factor,])
Enrichment[5,2] <- nrow(ExpressedRegulated)
Enrichment[6,1] <- length(Expressed[ Expressed$Factor %in% GO$Factor,])
Enrichment[6,2] <- nrow(Expressed)
Enrichment[7,1] <- length(TFs[ TFs[,1] %in% GO$Factor,])
Enrichment[7,2] <- nrow(TFs)

# Plot it
par(mfcol=c(1,1))
barplot(Enrichment[,1]/Enrichment[,2], las=2, ylim=c(0,0.2), ylab="Precision", names=c("IMAGE","dMotif+dExprs","Dynamic motifs","All motifs","Dynamic expression","All expressed","All TFs"), col=c("green","purple","red","pink","blue","lightblue","grey"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)