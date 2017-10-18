```R
# Source the necessary libaries
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(GenomicRanges)

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

## Setup to capture results
MeanDF <- data.frame(matrix(ncol=2, nrow=2))
SDDF <- data.frame(matrix(ncol=2, nrow=2))

## IMAGE
# Setup up to run it
Targets <- Target_Genes$PPARG_3.motif
NumberTargets <- nrow(Targets[ Targets$Target == 1,])

Genes <- Test
Genes$Top <- 0
Genes[ (Genes$Control_FDR <= 0.01 & Genes$Control_logFC > 0 & Genes$Control_logFC -1 >= (Genes$KD_logFC)) | (Genes$Control_FDR <= 0.01 & Genes$Control_logFC < 0 & (Genes$Control_logFC + 1) <= (Genes$KD_logFC)) ,"Top"] <- 1
Genes <- Genes[ order(Genes$Symbol_Human, -Genes$Top),]
Genes <- Genes[ duplicated(Genes$Symbol_Human)==F,]

# Calculate the precision
Precision <- nrow(Targets[ Targets$Target == 1 & Targets$Factor %in% Genes[ Genes$Top == 1,"Symbol_Human"],])/nrow(Genes[ Genes$Top == 1,])

# Randomize the targets, find overlap with significant genes
RandomTargets <- data.frame(matrix(ncol=1, nrow=1000))
for (m in 1:1000) {
	Targets$Random <- 0
	Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
	RandomTargets[m,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Genes[ Genes$Top == 1,"Symbol_Human"],])/nrow(Genes[ Genes$Top == 1,])
	}

# Randomize both groups for IMAGE
Randomized <- data.frame(matrix(ncol=1, nrow=1000))
for (q in 1:1000) {
	Genes$Random <- 0
	Genes[ sample(1:nrow(Genes),nrow(Genes[ Genes$Top == 1,])),"Random"] <- 1
	Targets$Random <- 0
	Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
	Randomized[q,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% Genes[ Genes$Random == 1,"Symbol_Human"],])/nrow(Genes[ Genes$Random == 1,])
	}

# Save the results	
MeanDF[1,] <- c(mean(Precision/Randomized[,1]), mean(mean(RandomTargets[,1])/Randomized[,1]))
SDDF[1,] <- c(sd(Precision/Randomized[,1]), sd(mean(RandomTargets[,1])/Randomized[,1]))

## PPARG ChIP
# Import ChIP-seq data
PPAR <- read.table("Data/PPARG.pos", quote="\"", stringsAsFactors=FALSE)
PeakGRange <- GRanges(seqnames = PPAR$V2, IRanges(start = PPAR$V3, end = PPAR$V4), strand = PPAR$V5)

# Import TSS and convert
TSS <- read.delim("Data/mm9.tss", header=F)
Convert <- read.delim("Data/Mouse_Human.txt", stringsAsFactors=FALSE)
colnames(Convert) <- c("Ensembl_Mouse","Ensembl_Human")
Conv1 <- as.data.frame(org.Mm.egREFSEQ)
Conv2 <- as.data.frame(org.Mm.egENSEMBL)
TSS <- merge(TSS, Conv1, by.x="V1", by.y="accession")
TSS <- merge(TSS, Conv2, by="gene_id")
TSS <- merge(TSS, Convert, by.x="ensembl_id", by.y="Ensembl_Mouse")
Conv1 <- as.data.frame(org.Hs.egENSEMBL)
Conv2 <- as.data.frame(org.Hs.egSYMBOL)
Conv <- merge(Conv1, Conv2, by="gene_id")
colnames(Conv) <- c("GeneID","Ensembl_Human","Symbol_Human")
TSS <- merge(TSS, Conv, by="Ensembl_Human")

# Define target genes
Size <- 25000
TSS$V3 <- TSS$V3 - Size
TSS$V4 <- TSS$V4 + Size
GeneGRange <- GRanges(seqnames = TSS$V2, IRanges(start = TSS$V3, end = TSS$V4), strand = rep("+", nrow(TSS)))
overlap <- as.data.frame(findOverlaps(GeneGRange, PeakGRange))
overlap <- overlap[ duplicated(overlap[,1]) == F,]
ChIPTargets <- TSS[ overlap$queryHits,] 
ChIPTargets <- ChIPTargets[ duplicated(ChIPTargets$Symbol_Human)==F,]

# Find the precision
Precision <- nrow(Targets[ Targets$Target == 1 & Targets$Factor %in% ChIPTargets$Symbol_Human,])/nrow(ChIPTargets)

# Randomize the targets, find overlap with enhancers
RandomTargets <- data.frame(matrix(ncol=1, nrow=1000))
for (m in 1:1000) {
	Targets$Random <- 0
	Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
	RandomTargets[m,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% ChIPTargets$Symbol_Human,])/nrow(ChIPTargets)
	}

# Randomize both groups for IMAGE
Randomized <- data.frame(matrix(ncol=1, nrow=1000))
for (q in 1:1000) {
	Genes$Random <- 0
	Genes[ sample(1:nrow(Genes),10000),"Random"] <- 1
	Targets$Random <- 0
	Targets[ sample(1:nrow(Targets), NumberTargets),"Random"] <- 1
	Randomized[q,1] <- nrow(Targets[ Targets$Random == 1 & Targets$Factor %in% ChIPTargets$Symbol_Human,])/nrow(Genes[ Genes$Random == 1,])
	}
	
# Save the results
MeanDF[2,] <- c(mean(Precision/Randomized[,1]), mean(mean(RandomTargets[,1])/Randomized[,1]))
SDDF[2,] <- c(sd(Precision/Randomized[,1]), sd(mean(RandomTargets[,1])/Randomized[,1]))

## Plot the results
# Setup
MeanDF <- t(as.matrix(MeanDF))
SDDF <- t(as.matrix(SDDF))

# Plot it
par(mfcol=c(1,1))
B <- barplot(log2(MeanDF), beside=T, names=c("IMAGE","Shuffled","ChIP","Shuffled"), las=2, ylab="log2 Enrichment", ylim=c(-0.1,2.5), col=c("green","grey","blue","grey"))
arrows(B, log2(MeanDF+SDDF), B, log2(MeanDF-SDDF), code=3, angle=90, length=0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)