```R
# Source the necessary libraries
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Import the data
load("Data/3T3.R")

# Import functions for RSA analysis
source("Data/OPI_Functions.R")

## Process raw data
# Import the lists of transcription factors
TFs <- read.table("Data/TFs.txt", quote="\"", comment.char="")
MotifToFactor <- read.delim("Data/Genename_Motif.txt", header=FALSE)

# Collaps hits
Hits <- Result[ Result$CausalTF == 1,]
Hits <- Hits[ duplicated(Hits$Factor)==F,]

## Import the ISMARA results (All with Z >= 2)
MARA <- read.delim("Data/ISMARA", header=F)
MARA[,1] <- toupper(MARA[,1])
MARA <- MARA[ MARA$V1 %in% TFs$V1,]
MARA <- MARA[ duplicated(MARA[,1]) == F,]

## Make a data frame to save the results
Impact <- data.frame(matrix(ncol=2, nrow=15))

## Process GO data
# Import and process data
GO <- read.delim("Data/FatCell_GO.txt", header=TRUE, stringsAsFactors=FALSE)

# Save enrichments
Impact[1,1] <- nrow(GO[ GO$Factor %in% Hits[ Hits$Pearsons >= 0.8 ,"Factor"],])
Impact[1,2] <- nrow(Hits[ Hits$Pearsons >= 0.8 ,])
Impact[2,1] <- nrow(GO[ GO$Factor %in% Hits$Factor,])
Impact[2,2] <- nrow(Hits)
Impact[3,1] <- nrow(GO[ GO$Factor %in% MARA[ MARA$V2 >= 0.8,"V1"],])
Impact[3,2] <- nrow(MARA[ MARA$V2 >= 0.8,])
Impact[4,1] <- nrow(GO[ GO$Factor %in% MARA$V1,])
Impact[4,2] <- nrow(MARA)
Impact[5,1] <- nrow(GO[ GO$Factor %in% TFs$V1,])
Impact[5,2] <- nrow(TFs)

## Overexpression screen
# Import and process data
Convert <- read.delim("Data/Mouse_Human.txt", stringsAsFactors=FALSE)
Screen <- read.delim("Data/Screen.txt", header=TRUE, stringsAsFactors=FALSE)
Factors <- read.table("Data/ScreenedFactors.txt", quote="\"", comment.char="")
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

# Save enrichments
Impact[6,1] <- nrow(Screen[ Screen$Symbol_Human %in% Hits[ Hits$Pearsons >= 0.8, "Factor"],])
Impact[6,2] <- nrow(Factors[ Factors$Symbol_Human %in% Hits[ Hits$Pearsons >= 0.8, "Factor"],])
Impact[7,1] <- nrow(Screen[ Screen$Symbol_Human %in% Hits[ , "Factor"],])
Impact[7,2] <- nrow(Factors[ Factors$Symbol_Human %in% Hits[ , "Factor"],])
Impact[8,1] <- nrow(Screen[ Screen$Symbol_Human %in% MARA[ MARA$V2 >= 0.8, "V1"],])
Impact[8,2] <- nrow(Factors[ Factors$Symbol_Human %in% MARA[ MARA$V2 >= 0.8, "V1"],])
Impact[9,1] <- nrow(Screen[ Screen$Symbol_Human %in% MARA[ , "V1"],])
Impact[9,2] <- nrow(Factors[ Factors$Symbol_Human %in% MARA[ , "V1"],])
Impact[10,1] <- nrow(Screen[ Screen$Symbol_Human %in% TFs$V1,])
Impact[10,2] <- nrow(Factors[ Factors$Symbol_Human %in% TFs$V1,])

## Knockdown screen
# Import and process data
Screen <- read.delim("Data/WinnefeldScreen.txt")
Screen <- Screen[ Screen$GeneSymbol %in% TFs$V1,]
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

# Save enrichments
Impact[11,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% Hits[ Hits$Pearsons >= 0.8, "Factor"] & ResultScreen$Hit == 1,])
Impact[11,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% Hits[ Hits$Pearsons >= 0.8, "Factor"],])
Impact[12,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% Hits[ , "Factor"] & ResultScreen$Hit == 1,])
Impact[12,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% Hits[ , "Factor"],])
Impact[13,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% MARA[ MARA$V2 >= 0.8, "V1"] & ResultScreen$Hit == 1,])
Impact[13,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% MARA[ MARA$V2 >= 0.8, "V1"],])
Impact[14,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% MARA[ , "V1"] & ResultScreen$Hit == 1,])
Impact[14,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% MARA[ , "V1"],])
Impact[15,1] <- nrow(ResultScreen[ ResultScreen$Factor %in% TFs$V1 & ResultScreen$Hit == 1,])
Impact[15,2] <- nrow(ResultScreen[ ResultScreen$Factor %in% TFs$V1,])

## Process all enrichments and plot it
# Calculate enrichments
Impact[,3] <- Impact[,1] / Impact[,2]
for (i in 1:4) { Impact[i,3] <- Impact[i,3] / Impact[5,3] }
for (i in 6:9) { Impact[i,3] <- Impact[i,3] / Impact[10,3] }
for (i in 11:14) { Impact[i,3] <- Impact[i,3] / Impact[15,3] }
Impact[,3] <- log2(Impact[,3])

# Make the plots
par(mfcol=c(1,3))
barplot(Impact[c(1,3,2,4),3], las=2, ylab="Log2 enrichment", names=c("IMAGE","MARA","IMAGE","MARA"), main="Gene ontology", col=c("green","orange", "green","orange"))
mtext("Correlating", at = 1.25, side=1, line=2.5)
mtext("All", at = 3.75, side=1, line=2.5)

barplot(Impact[c(6,8,7,9),3], las=2, ylab="Log2 enrichment", names=c("IMAGE","MARA","IMAGE","MARA"), main="Overexpression screen", col=c("green","orange", "green","orange"))
mtext("Correlating", at = 1.25, side=1, line=2.5)
mtext("All", at = 3.75, side=1, line=2.5)

barplot(Impact[c(11,13,12,14),3], las=2, ylab="Log2 enrichment", names=c("IMAGE","MARA","IMAGE","MARA"), main="Knockdown screen", col=c("green","orange", "green","orange"))
mtext("Correlating", at = 1.25, side=1, line=2.5)
mtext("All", at = 3.75, side=1, line=2.5)
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)