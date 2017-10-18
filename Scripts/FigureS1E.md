```R
# Source the necessary libaries
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Import the data
load("Data/3T3.R")

# Import the screening results, conversion files and transcription factors
Screen <- read.delim("Data/Screen.txt", header=TRUE, stringsAsFactors=FALSE)
Factors <- read.table("Data/ScreenedFactors.txt", quote="\"", comment.char="")
TFs <- read.table("Data/TFs.txt", quote="\"", comment.char="")
Convert <- read.delim("Data/Mouse_Human.txt", stringsAsFactors=FALSE)
MotifToFactor <- read.delim("Data/Genename_Motif.txt", header=FALSE)

# Convert symbols between mouse and man of the screening results
colnames(Convert) <- c("Ensembl_Mouse","Ensembl_Human")
Conv1 <- as.data.frame(org.Hs.egENSEMBL)
Conv2 <- as.data.frame(org.Hs.egSYMBOL)
Conv <- merge(Conv1, Conv2, by="gene_id")
colnames(Conv) <- c("GeneID","Ensembl_Human","Symbol_Human")
Convert <- merge(Convert, Conv[,c(2,3)], by="Ensembl_Human")
Conv1 <- as.data.frame(org.Mm.egENSEMBL)
Conv2 <- as.data.frame(org.Mm.egSYMBOL)
Conv <- merge(Conv1, Conv2, by="gene_id")
colnames(Conv) <- c("GeneID","Ensembl_Mouse","Symbol_Mouse")
Convert <- merge(Convert, Conv[,c(2,3)], by="Ensembl_Mouse")
Screen <- merge(Screen, Convert[,c(3,4)], by.x="Gene", by.y="Symbol_Mouse")
Screen <- Screen[ duplicated(Screen$Symbol_Human)==F,]
Screen <- Screen[ (Screen$Pvalue1 <= 0.05 & Screen$Pvalue2 <= 0.05) | (Screen$Pvalue1 <= 0.05 & Screen$Pvalue3 <= 0.05)| (Screen$Pvalue3 <= 0.05 & Screen$Pvalue2 <= 0.05), ]
Factors <- merge(Factors, Convert[,c(3,4)], by.x="V1", by.y="Symbol_Mouse")
Factors <- Factors[ duplicated(Factors$Symbol_Human)==F,]
rm(Conv1)
rm(Conv2)
rm(Conv)
rm(Convert)

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

# Calculate mean, SEM and p-value for each group
Impact <- data.frame(matrix(ncol=2, nrow=8))
Impact[1,1] <- nrow(Screen[ Screen$Symbol_Human %in% Result[ Result$CausalTF == 1, "Factor"],])
Impact[1,2] <- nrow(Factors[ Factors$Symbol_Human %in% Result[ Result$CausalTF == 1, "Factor"],])
Impact[2,1] <- nrow(Screen[ Screen$Symbol_Human %in% Combined$V1,])
Impact[2,2] <- nrow(Factors[ Factors$Symbol_Human %in% Combined$V1,])
Impact[3,1] <- nrow(Screen[ Screen$Symbol_Human %in% MotifsRegulated$V1,])
Impact[3,2] <- nrow(Factors[ Factors$Symbol_Human %in% MotifsRegulated$V1,])
Impact[4,1] <- nrow(Screen[ Screen$Symbol_Human %in% MotifsAll$V1,])
Impact[4,2] <- nrow(Factors[ Factors$Symbol_Human %in% MotifsAll$V1,])
Impact[5,1] <- nrow(Screen[ Screen$Symbol_Human %in% ExpressedRegulated$Factor,])
Impact[5,2] <- nrow(Factors[ Factors$Symbol_Human %in% ExpressedRegulated$Factor,])
Impact[6,1] <- nrow(Screen[ Screen$Symbol_Human %in% Expressed$Factor,])
Impact[6,2] <- nrow(Factors[ Factors$Symbol_Human %in% Expressed$Factor,])
Impact[8,1] <- nrow(Screen[ Screen$Symbol_Human %in% TFs$V1,])
Impact[8,2] <- nrow(Factors[ Factors$Symbol_Human %in% TFs$V1,])

# Perform permutations on random selection
RandomResult <- data.frame(matrix(ncol=2, nrow=1000))
for (i in 1:1000) {
  Random <- data.frame(Factor = TFs[ sample(1:nrow(TFs),500),1])
  RandomResult[i,1] <- nrow(Screen[ Screen$Symbol_Human %in% Random$Factor,])
  RandomResult[i,2] <- nrow(Factors[ Factors$Symbol_Human %in% Random$Factor,])
}

Impact[7,1] <- mean(RandomResult[,1])
Impact[7,2] <- mean(RandomResult[,2])

# Calculate enrichments
Impact[,3] <- Impact[,1] / Impact[,2]
Impact[,3] <- Impact[,3] / Impact[8,3]
Impact[,3] <- log2(Impact[,3])

# Make the plot
barplot(Impact[c(7:1),3], las=2, ylim=c(0,0.6), ylab="Log2 Enrichment over all TFs", names=c("Random","All","Dynamic","All","Dynamic","dMotif+cExprs","IMAGE"), col=c("grey","blue","blue","red","red","purple","green"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)