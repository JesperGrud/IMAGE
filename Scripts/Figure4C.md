```R
# Import functions for RSA analysis
source("Data/OPI_Functions.R")

# Import the data
load("Data/3T3.R")

## Define groups of interest
# Import the list of transcription factors
TFs <- read.table("Data/TFs.txt", quote="\"", comment.char="")
MotifToFactor <- read.delim("Data/Genename_Motif.txt", header=FALSE)

## Make groups of transcription factors based on expression patterns
# All expressed
Expressed <- data.frame(Factor = TFs[ TFs[,1] %in% GeneRPKM[ GeneRPKM$Max >= 0,"Factor"],])
Expressed <- data.frame(Factor = Expressed[ duplicated(Expressed$Factor)==F,])

# Regulated during differentiation
ExpressedRegulated <- data.frame(Factor = TFs[ TFs[,1] %in% GeneRPKM[ GeneRPKM$Max >= 0 & ((GeneRPKM$`FDR_Cond3-vs-Cond2` <= 0.05) | (GeneRPKM$`FDR_Cond4-vs-Cond2` <= 0.05) | (GeneRPKM$`FDR_Cond4-vs-Cond3` <= 0.05) | (GeneRPKM$`FDR_Cond4-vs-Cond1` <= 0.05)  | (GeneRPKM$`FDR_Cond3-vs-Cond1` <= 0.05)  | (GeneRPKM$`FDR_Cond2-vs-Cond1` <= 0.05)),"Factor"],])

# Make groups of transcription factors based on motif enrichments
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
rm(MotifsInd)
rm(MotifsRep)

# Make groups based on combinations of motifs and expression
Combined <- MotifsRegulated[ MotifsRegulated$V1 %in% ExpressedRegulated$Factor,]


## Import the screen data
Screen <- read.delim("Data/WinnefeldScreen.txt")
Screen <- Screen[ Screen$GeneSymbol %in% TFs$V1,]

## Setup and run RSA
opts = list(LB=0.8,UB=1.96,outputFile="RSA_out",inputFile=NA,reverse=TRUE,bonferroni=1);
t = data.frame(Gene_ID = Screen$GeneSymbol, Well_ID = Screen$well, Score = Screen$score)
reverse = OPI(t$Gene_ID,t$Score,opts,t)
reverse$Induce <- 0
reverse[ reverse$OPI_Hit == 1 & reverse$`#hitWell` >= 2, "Induce"] <- 1
reverse <- reverse[ duplicated(reverse$Gene_ID)==F,]
t = data.frame(Gene_ID = Screen$GeneSymbol, Well_ID = Screen$well, Score = -1*Screen$score)
forward = OPI(t$Gene_ID,t$Score,opts,t)
forward$Repress <- 0
forward[ forward$OPI_Hit == 1 & forward$`#hitWell` >= 2, "Repress"] <- 1
forward <- forward[ duplicated(forward$Gene_ID)==F,]
ResultScreen <- merge(forward[,c("Gene_ID","Repress")], reverse[,c("Gene_ID","Induce")], by="Gene_ID")
colnames(ResultScreen)[1] <- "Factor"
ResultScreen$Hit <- 0
ResultScreen[ ResultScreen$Repress == 1 | ResultScreen$Induce == 1, "Hit"] <- 1

# Calculate the results
Impact <- data.frame(matrix(ncol=2, nrow=8))
Impact[1,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% Result[ Result$CausalTF == 1, "Factor"] & ResultScreen$Hit == 1,])
Impact[1,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% Result[ Result$CausalTF == 1, "Factor"],])
Impact[2,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% Combined$V1 & ResultScreen$Hit == 1,])
Impact[2,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% Combined$V1,])
Impact[3,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% MotifsRegulated$V1 & ResultScreen$Hit == 1,])
Impact[3,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% MotifsRegulated$V1,])
Impact[4,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% MotifsAll$V1 & ResultScreen$Hit == 1,])
Impact[4,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% MotifsAll$V1,])
Impact[5,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% ExpressedRegulated$Factor & ResultScreen$Hit == 1,])
Impact[5,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% ExpressedRegulated$Factor,])
Impact[6,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% Expressed$Factor & ResultScreen$Hit == 1,])
Impact[6,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% Expressed$Factor,])
Impact[8,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% TFs$V1 & ResultScreen$Hit == 1,])
Impact[8,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% TFs$V1,])

# Perform permutations on random selection
RandomResult <- data.frame(matrix(ncol=2, nrow=1000))
for (i in 1:1000) {
Random <- data.frame(Factor = TFs[ sample(1:nrow(TFs),500),1])
RandomResult[i,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% Random$Factor & ResultScreen$Hit == 1,])
RandomResult[i,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% Random$Factor,])
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
[Back to overview of Figure 4](../Links/Figure4.md)