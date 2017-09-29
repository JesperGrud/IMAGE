```R
# Import libraries
library(edgeR)
library(pheatmap)

# Import the data
Counts <- read.delim("Data/Candidates.count.txt")

# Analyze using edgeR
rownames(Counts) <- Counts$RefSeq
DGE <- DGEList(Counts[,c(9:50)], group = substr(colnames(Counts)[c(9:50)],0,regexpr("SM",colnames(Counts)[c(9:50)])-5))
DGE <- calcNormFactors(DGE)
DGE <- estimateCommonDisp(DGE)
DGE <- estimateTagwiseDisp(DGE)

# Get RPKM
Normalized <- rpkm(DGE, gene.length = Counts$countLength, normalized.lib.sizes = T, log = T)
Normalized <- as.data.frame(Normalized)

# Calculate average for each replicate
for (i in 1:42) { Normalized[,i] <- 2^Normalized[,i] }
Normalized$HSF1_D0 <- rowMeans(Normalized[,c(1,2)])
Normalized$HSF1_D1 <- rowMeans(Normalized[,c(3,4)])
Normalized$HSF1_D7 <- rowMeans(Normalized[,c(5,6)])
Normalized$MAZ_D0 <- rowMeans(Normalized[,c(7,8)])
Normalized$MAZ_D1 <- rowMeans(Normalized[,c(9,10)])
Normalized$MAZ_D7 <- rowMeans(Normalized[,c(11,12)])
Normalized$MYBL1_D0 <- rowMeans(Normalized[,c(13,14)])
Normalized$MYBL1_D1 <- rowMeans(Normalized[,c(15,16)])
Normalized$MYBL1_D7 <- rowMeans(Normalized[,c(17,18)])
Normalized$NFIL3_D0 <- rowMeans(Normalized[,c(19,20)])
Normalized$NFIL3_D1 <- rowMeans(Normalized[,c(21,22)])
Normalized$NFIL3_D7 <- rowMeans(Normalized[,c(23,24)])
Normalized$SATB1_D0 <- rowMeans(Normalized[,c(25,26)])
Normalized$SATB1_D1 <- rowMeans(Normalized[,c(27,28)])
Normalized$SATB1_D7 <- rowMeans(Normalized[,c(29,30)])
Normalized$TSHZ1_D0 <- rowMeans(Normalized[,c(31,32)])
Normalized$TSHZ1_D1 <- rowMeans(Normalized[,c(33,34)])
Normalized$TSHZ1_D7 <- rowMeans(Normalized[,c(35,36)])
Normalized$UNI_D0 <- rowMeans(Normalized[,c(37,38)])
Normalized$UNI_D1 <- rowMeans(Normalized[,c(39,40)])
Normalized$UNI_D7 <- rowMeans(Normalized[,c(41,42)])
Normalized$RefSeq <- rownames(Normalized)

# Compare within samples over time
UNI_D0_D1 <- as.data.frame(topTags(exactTest(DGE, pair = c("UNI_D0","UNI_D1")), n = nrow(Counts)))
UNI_D0_D1$RefSeq <- rownames(UNI_D0_D1)
UNI_D0_D7 <- as.data.frame(topTags(exactTest(DGE, pair = c("UNI_D0","UNI_D7")), n = nrow(Counts)))
UNI_D0_D7$RefSeq <- rownames(UNI_D0_D7)

HSF_D0_D1 <- as.data.frame(topTags(exactTest(DGE, pair = c("HSF1_D0","HSF1_D1")), n = nrow(Counts)))
HSF_D0_D1$RefSeq <- rownames(HSF_D0_D1)
HSF_D0_D7 <- as.data.frame(topTags(exactTest(DGE, pair = c("HSF1_D0","HSF1_D7")), n = nrow(Counts)))
HSF_D0_D7$RefSeq <- rownames(HSF_D0_D7)

MAZ_D0_D1 <- as.data.frame(topTags(exactTest(DGE, pair = c("MAZ_D0","MAZ_D1")), n = nrow(Counts)))
MAZ_D0_D1$RefSeq <- rownames(MAZ_D0_D1)
MAZ_D0_D7 <- as.data.frame(topTags(exactTest(DGE, pair = c("MAZ_D0","MAZ_D7")), n = nrow(Counts)))
MAZ_D0_D7$RefSeq <- rownames(MAZ_D0_D7)

MYBL1_D0_D1 <- as.data.frame(topTags(exactTest(DGE, pair = c("MYBL1_D0","MYBL1_D1")), n = nrow(Counts)))
MYBL1_D0_D1$RefSeq <- rownames(MYBL1_D0_D1)
MYBL1_D0_D7 <- as.data.frame(topTags(exactTest(DGE, pair = c("MYBL1_D0","MYBL1_D7")), n = nrow(Counts)))
MYBL1_D0_D7$RefSeq <- rownames(MYBL1_D0_D7)

NFIL3_D0_D1 <- as.data.frame(topTags(exactTest(DGE, pair = c("NFIL3_D0","NFIL3_D1")), n = nrow(Counts)))
NFIL3_D0_D1$RefSeq <- rownames(NFIL3_D0_D1)
NFIL3_D0_D7 <- as.data.frame(topTags(exactTest(DGE, pair = c("NFIL3_D0","NFIL3_D7")), n = nrow(Counts)))
NFIL3_D0_D7$RefSeq <- rownames(NFIL3_D0_D7)

SATB1_D0_D1 <- as.data.frame(topTags(exactTest(DGE, pair = c("SATB1_D0","SATB1_D1")), n = nrow(Counts)))
SATB1_D0_D1$RefSeq <- rownames(SATB1_D0_D1)
SATB1_D0_D7 <- as.data.frame(topTags(exactTest(DGE, pair = c("SATB1_D0","SATB1_D7")), n = nrow(Counts)))
SATB1_D0_D7$RefSeq <- rownames(SATB1_D0_D7)

TSHZ1_D0_D1 <- as.data.frame(topTags(exactTest(DGE, pair = c("TSHZ1_D0","TSHZ1_D1")), n = nrow(Counts)))
TSHZ1_D0_D1$RefSeq <- rownames(TSHZ1_D0_D1)
TSHZ1_D0_D7 <- as.data.frame(topTags(exactTest(DGE, pair = c("TSHZ1_D0","TSHZ1_D7")), n = nrow(Counts)))
TSHZ1_D0_D7$RefSeq <- rownames(TSHZ1_D0_D7)

# Combine results
Control <- merge(UNI_D0_D1[,c(1,4,5)], UNI_D0_D7[,c(1,4,5)], by="RefSeq")
HSF <- merge(HSF_D0_D1[,c(1,4,5)], HSF_D0_D7[,c(1,4,5)], by="RefSeq")
MYBL1 <- merge(MYBL1_D0_D1[,c(1,4,5)], MYBL1_D0_D7[,c(1,4,5)], by="RefSeq")
NFIL3 <- merge(NFIL3_D0_D1[,c(1,4,5)], NFIL3_D0_D7[,c(1,4,5)], by="RefSeq")
SATB1 <- merge(SATB1_D0_D1[,c(1,4,5)], SATB1_D0_D7[,c(1,4,5)], by="RefSeq")
MAZ <- merge(MAZ_D0_D1[,c(1,4,5)], MAZ_D0_D7[,c(1,4,5)], by="RefSeq")
TSHZ1 <- merge(TSHZ1_D0_D1[,c(1,4,5)], TSHZ1_D0_D7[,c(1,4,5)], by="RefSeq")

# Plot the expression of adipocyte marker genes
Markers <- c("LPL","ADIPOQ","PPARG","PLIN1","PCK1","FABP4","LIPE", "SCD1","CIDEC","CEBPA", "PLIN4")
HeatmapData <- merge(Normalized, Counts[,c("RefSeq","Symbol")], by="RefSeq")
HeatmapData <- HeatmapData[ HeatmapData$Symbol %in% Markers,]
rownames(HeatmapData) <- HeatmapData$Symbol
pheatmap((HeatmapData[,c(2:43)]), cluster_cols = F, scale="row")
```

[Back to start](../README.md)<br>
[Back to overview of Figure 5](../Links/Figure5.md)