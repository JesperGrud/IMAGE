```R
# Source the necessary libaries
library(org.Mm.eg.db)
library(org.Hs.eg.db)

# Import the data
load("Data/3T3.R")

# Import the knockdown data
Test <- read.delim("Data/Adipogenesis_GSE14004.txt")
Test2 <- read.delim("Data/Adipogenesis_siPPARg_GSE14004.txt")
colnames(Test)[2] <- "Control_FDR"
colnames(Test)[6] <- "Control_logFC"
colnames(Test)[7] <- "Symbol"
colnames(Test2)[7] <- "Symbol"
colnames(Test2)[6] <- "KD_logFC"
colnames(Test2)[2] <- "KD_FDR"
Test <- merge(Test[,c(1,2,6)], Test2[,c(1,2,6,7)], by="ID")

# Convert to human IDs
Convert <- read.delim("Data/Mouse_Human.txt", stringsAsFactors=FALSE)
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
Test <- merge(Test, Convert[,c(3,4)], by.x="Symbol", by.y="Symbol_Mouse")
Test <- Test[ duplicated(Test$Symbol_Human)==F,]
rm(Conv1)
rm(Conv2)
rm(Conv)
rm(Convert)

# Import ChIP-seq data and define target genes
PPAR <- read.table("Data/PPARG.pos", quote="\"", stringsAsFactors=FALSE)
PPAR <- PPAR[ PPAR$V11 >= 10,]
PeakGRange <- GRanges(seqnames = PPAR$V2, IRanges(start = PPAR$V3, end = PPAR$V4), strand = PPAR$V5)
GeneRPKM$TSS <- 0
for (i in 1:nrow(GeneRPKM)) { if (GeneRPKM[i,"Strand"] == "+") { GeneRPKM[i,"TSS"] <- GeneRPKM[i,"Start"] } else { GeneRPKM[i,"TSS"] <- GeneRPKM[i,"End"] }  }
GeneRPKM$RegionStart <- GeneRPKM$TSS - 10000
GeneRPKM$RegionEnd <- GeneRPKM$TSS + 10000
GeneRPKM[ GeneRPKM$RegionStart < 0 , "RegionStart"] <- 0
GeneGRange <- GRanges(seqnames = GeneRPKM$Chr, IRanges(start = GeneRPKM$RegionStart, end = GeneRPKM$RegionEnd), strand = rep("+", nrow(GeneRPKM)))
overlap <- as.data.frame(findOverlaps(GeneGRange, PeakGRange))
overlap <- overlap[ duplicated(overlap[,1]) == F,]
ChIPTargets <- GeneRPKM[ overlap$queryHits,]

# Calculate enrichment of KD affected genes
ResultEnrichment <- data.frame(matrix(ncol=2, nrow=6))
Test <- Test[ Test$Control_FDR <= 0.01,]
Affected <- Test[ (Test$Control_logFC > 0 & Test$Control_logFC -1 >= (Test$KD_logFC)) | (Test$Control_logFC < 0 & (Test$Control_logFC + 1) <= (Test$KD_logFC)) ,]
Affected <- Affected[ duplicated(Affected$Symbol_Human)==F,]
Prediction <- Target_Genes$PPARG_1.motif
Predicted <- Prediction[ Prediction$Target == 1, ]
Predicted <- Predicted[ duplicated(Predicted$Factor)==F,]
ResultEnrichment[1,1] <- nrow(Predicted[ Predicted$Factor %in% Affected$Symbol_Human,])
ResultEnrichment[1,2] <- nrow(Predicted[ Predicted$Factor %in% Test$Symbol_Human,])
ResultEnrichment[2,1] <- nrow(ChIPTargets[ ChIPTargets$Factor %in% Affected$Symbol_Human,])
ResultEnrichment[2,2] <- nrow(ChIPTargets[ ChIPTargets$Factor %in% Test$Symbol_Human,])
ResultEnrichment[6,1] <- nrow(Prediction[ Prediction$Factor %in% Affected$Symbol_Human,])
ResultEnrichment[6,2] <- nrow(Prediction[ Prediction$Factor %in% Test$Symbol_Human,])

# Perform 1000 permutations of random selection (both random genes, random TFs, and random enhancers)
Randomize <- data.frame(matrix(ncol=6, nrow=1000))
for (i in 1:1000) {
# Random prediction
Random <- Prediction[ sample(1:nrow(Prediction),nrow(Predicted)),]
Randomize[i,1] <- nrow(Random[ Random$Factor %in% Affected$Symbol_Human,])
Randomize[i,2] <- nrow(Random[ Random$Factor %in% Test$Symbol_Human,])
# Random TF
Random <- Target_Genes[[sample(1:length(Target_Genes),1)]]
Random <- Random[ Random$Target == 1,]
Randomize[i,3] <- nrow(Random[ Random$Factor %in% Affected$Symbol_Human,])
Randomize[i,4] <- nrow(Random[ Random$Factor %in% Test$Symbol_Human,])
# Random enhancers
RandoEnhancers <- Enhancers[ sample(1:nrow(Enhancers), nrow(PPAR)),]
PeakGRange <- GRanges(seqnames = RandoEnhancers$V1, IRanges(start = RandoEnhancers$V2, end = RandoEnhancers$V3), strand = rep("+",nrow(RandoEnhancers)))
overlap <- as.data.frame(findOverlaps(GeneGRange, PeakGRange))
overlap <- overlap[ duplicated(overlap[,1]) == F,]
RandoTargets <- GeneRPKM[ overlap$queryHits,]
Randomize[i,5] <- nrow(RandoTargets[ RandoTargets$Factor %in% Affected$Symbol_Human,])
Randomize[i,6] <- nrow(RandoTargets[ RandoTargets$Factor %in% Test$Symbol_Human,])
}

ResultEnrichment[3,1] <- mean(Randomize[,3])
ResultEnrichment[3,2] <- mean(Randomize[,4])
ResultEnrichment[4,1] <- mean(Randomize[,1])
ResultEnrichment[4,2] <- mean(Randomize[,2])
ResultEnrichment[5,1] <- mean(Randomize[,5])
ResultEnrichment[5,2] <- mean(Randomize[,6])

# Calculate enrichment and make the plot
par(mfcol=c(1,1))
ResultEnrichment[,3] <- ResultEnrichment[,1] / ResultEnrichment[,2]
ResultEnrichment[,3] <- ResultEnrichment[,3] / ResultEnrichment[6,3]
ResultEnrichment[,3] <- log2(ResultEnrichment[,3])
barplot(ResultEnrichment[c(2,5,1,3,4),3], ylim=c(-0.1,1), las=1, ylab="Enrichment", col=c("blue","grey","green","grey","grey"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)