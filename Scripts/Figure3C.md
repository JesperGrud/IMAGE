```R
# Source the necessary libaries
library(GenomicRanges)
library(pheatmap)

# Import the data
load("Data/3T3.R")

# Import the site counts and peak calls
Jun <- read.table("Data/cJun.pos", quote="\"", stringsAsFactors=FALSE)
Junb <- read.table("Data/Junb.pos", quote="\"", stringsAsFactors=FALSE)
PPAR <- read.table("/Data/PPARG.pos", quote="\"", stringsAsFactors=FALSE)
GR <- read.table("Data/GR.pos", quote="\"", stringsAsFactors=FALSE)
Counts <- read.delim("Data/Site.counts")

## Fractional overlap
# Make a data frame for capturing the results
FractionResult <- data.frame(matrix(ncol=5, nrow=4))
# PPARg
All <- Target_Sites$PPARG_1.motif
Motif <- All[ All$NMotif > 0,]
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
AntiModel <- All[ All$NMotif > 0 & All$Contribution < 0,]
# Make GRange objects
PeakGRange <- GRanges(seqnames = PPAR$V2, IRanges(start = PPAR$V3, end = PPAR$V4), strand = PPAR$V5)
AllGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% All$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% All$PeakID, 3], end = Counts[ Counts$PeakID %in% All$PeakID, 4]), strand = Counts[ Counts$PeakID %in% All$PeakID, 5])
MotifGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% Motif$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% Motif$PeakID, 3], end = Counts[ Counts$PeakID %in% Motif$PeakID, 4]), strand = Counts[ Counts$PeakID %in% Motif$PeakID, 5])
ModelGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% Model$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% Model$PeakID, 3], end = Counts[ Counts$PeakID %in% Model$PeakID, 4]), strand = Counts[ Counts$PeakID %in% Model$PeakID, 5])
AntimodelGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% AntiModel$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% AntiModel$PeakID, 3], end = Counts[ Counts$PeakID %in% AntiModel$PeakID, 4]), strand = Counts[ Counts$PeakID %in% AntiModel$PeakID, 5])
# Get overlap
overlapAll <- as.data.frame(findOverlaps(PeakGRange, AllGRange))
overlapMotif <- as.data.frame(findOverlaps(PeakGRange, MotifGRange))
overlapModel <- as.data.frame(findOverlaps(PeakGRange, ModelGRange))
overlapAntimodel <- as.data.frame(findOverlaps(PeakGRange, AntimodelGRange))
# Remove duplicates
overlapAll <- overlapAll[ duplicated(overlapAll[,2]) == F,]
overlapMotif <- overlapMotif[ duplicated(overlapMotif[,2]) == F,]
overlapModel <- overlapModel[ duplicated(overlapModel[,2]) == F,]
overlapAntimodel <- overlapAntimodel[ duplicated(overlapAntimodel[,2]) == F,]
# Get fractions
FractionResult[1,1] <-"PPARg"
FractionResult[1,2] <- nrow(overlapAll)/nrow(All)
FractionResult[1,3] <- nrow(overlapMotif)/nrow(Motif)
FractionResult[1,4] <- nrow(overlapModel)/nrow(Model)
FractionResult[1,5] <- nrow(overlapAntimodel)/nrow(AntiModel)
# Jun
All <- Target_Sites$JUN_1.motif
Motif <- All[ All$NMotif > 0,]
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
AntiModel <- All[ All$NMotif > 0 & All$Contribution < 0,]
# Make GRange objects
PeakGRange <- GRanges(seqnames = Jun$V2, IRanges(start = Jun$V3, end = Jun$V4), strand = Jun$V5)
AllGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% All$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% All$PeakID, 3], end = Counts[ Counts$PeakID %in% All$PeakID, 4]), strand = Counts[ Counts$PeakID %in% All$PeakID, 5])
MotifGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% Motif$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% Motif$PeakID, 3], end = Counts[ Counts$PeakID %in% Motif$PeakID, 4]), strand = Counts[ Counts$PeakID %in% Motif$PeakID, 5])
ModelGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% Model$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% Model$PeakID, 3], end = Counts[ Counts$PeakID %in% Model$PeakID, 4]), strand = Counts[ Counts$PeakID %in% Model$PeakID, 5])
AntimodelGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% AntiModel$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% AntiModel$PeakID, 3], end = Counts[ Counts$PeakID %in% AntiModel$PeakID, 4]), strand = Counts[ Counts$PeakID %in% AntiModel$PeakID, 5])
# Get overlap
overlapAll <- as.data.frame(findOverlaps(PeakGRange, AllGRange))
overlapMotif <- as.data.frame(findOverlaps(PeakGRange, MotifGRange))
overlapModel <- as.data.frame(findOverlaps(PeakGRange, ModelGRange))
overlapAntimodel <- as.data.frame(findOverlaps(PeakGRange, AntimodelGRange))
# Remove duplicates
overlapAll <- overlapAll[ duplicated(overlapAll[,2]) == F,]
overlapMotif <- overlapMotif[ duplicated(overlapMotif[,2]) == F,]
overlapModel <- overlapModel[ duplicated(overlapModel[,2]) == F,]
overlapAntimodel <- overlapAntimodel[ duplicated(overlapAntimodel[,2]) == F,]
# Get fractions
FractionResult[2,1] <-"Jun"
FractionResult[2,2] <- nrow(overlapAll)/nrow(All)
FractionResult[2,3] <- nrow(overlapMotif)/nrow(Motif)
FractionResult[2,4] <- nrow(overlapModel)/nrow(Model)
FractionResult[2,5] <- nrow(overlapAntimodel)/nrow(AntiModel)
# GR
All <- Target_Sites$JUNB_1.motif
Motif <- All[ All$NMotif > 0,]
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
AntiModel <- All[ All$NMotif > 0 & All$Contribution < 0,]
# Make GRange objects
PeakGRange <- GRanges(seqnames = Junb$V2, IRanges(start = Junb$V3, end = Junb$V4), strand = Junb$V5)
AllGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% All$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% All$PeakID, 3], end = Counts[ Counts$PeakID %in% All$PeakID, 4]), strand = Counts[ Counts$PeakID %in% All$PeakID, 5])
MotifGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% Motif$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% Motif$PeakID, 3], end = Counts[ Counts$PeakID %in% Motif$PeakID, 4]), strand = Counts[ Counts$PeakID %in% Motif$PeakID, 5])
ModelGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% Model$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% Model$PeakID, 3], end = Counts[ Counts$PeakID %in% Model$PeakID, 4]), strand = Counts[ Counts$PeakID %in% Model$PeakID, 5])
AntimodelGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% AntiModel$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% AntiModel$PeakID, 3], end = Counts[ Counts$PeakID %in% AntiModel$PeakID, 4]), strand = Counts[ Counts$PeakID %in% AntiModel$PeakID, 5])
# Get overlap
overlapAll <- as.data.frame(findOverlaps(PeakGRange, AllGRange))
overlapMotif <- as.data.frame(findOverlaps(PeakGRange, MotifGRange))
overlapModel <- as.data.frame(findOverlaps(PeakGRange, ModelGRange))
overlapAntimodel <- as.data.frame(findOverlaps(PeakGRange, AntimodelGRange))
# Remove duplicates
overlapAll <- overlapAll[ duplicated(overlapAll[,2]) == F,]
overlapMotif <- overlapMotif[ duplicated(overlapMotif[,2]) == F,]
overlapModel <- overlapModel[ duplicated(overlapModel[,2]) == F,]
overlapAntimodel <- overlapAntimodel[ duplicated(overlapAntimodel[,2]) == F,]
# Get fractions
FractionResult[3,1] <-"Junb"
FractionResult[3,2] <- nrow(overlapAll)/nrow(All)
FractionResult[3,3] <- nrow(overlapMotif)/nrow(Motif)
FractionResult[3,4] <- nrow(overlapModel)/nrow(Model)
FractionResult[3,5] <- nrow(overlapAntimodel)/nrow(AntiModel)
# GR
All <- Target_Sites$NR3C1_1.motif
Motif <- All[ All$NMotif > 0,]
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
AntiModel <- All[ All$NMotif > 0 & All$Contribution < 0,]
# Make GRange objects
PeakGRange <- GRanges(seqnames = GR$V2, IRanges(start = GR$V3, end = GR$V4), strand = GR$V5)
AllGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% All$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% All$PeakID, 3], end = Counts[ Counts$PeakID %in% All$PeakID, 4]), strand = Counts[ Counts$PeakID %in% All$PeakID, 5])
MotifGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% Motif$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% Motif$PeakID, 3], end = Counts[ Counts$PeakID %in% Motif$PeakID, 4]), strand = Counts[ Counts$PeakID %in% Motif$PeakID, 5])
ModelGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% Model$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% Model$PeakID, 3], end = Counts[ Counts$PeakID %in% Model$PeakID, 4]), strand = Counts[ Counts$PeakID %in% Model$PeakID, 5])
AntimodelGRange <- GRanges(seqnames = Counts[ Counts$PeakID %in% AntiModel$PeakID, 2], IRanges(start = Counts[ Counts$PeakID %in% AntiModel$PeakID, 3], end = Counts[ Counts$PeakID %in% AntiModel$PeakID, 4]), strand = Counts[ Counts$PeakID %in% AntiModel$PeakID, 5])
# Get overlap
overlapAll <- as.data.frame(findOverlaps(PeakGRange, AllGRange))
overlapMotif <- as.data.frame(findOverlaps(PeakGRange, MotifGRange))
overlapModel <- as.data.frame(findOverlaps(PeakGRange, ModelGRange))
overlapAntimodel <- as.data.frame(findOverlaps(PeakGRange, AntimodelGRange))
# Remove duplicates
overlapAll <- overlapAll[ duplicated(overlapAll[,2]) == F,]
overlapMotif <- overlapMotif[ duplicated(overlapMotif[,2]) == F,]
overlapModel <- overlapModel[ duplicated(overlapModel[,2]) == F,]
overlapAntimodel <- overlapAntimodel[ duplicated(overlapAntimodel[,2]) == F,]
# Get fractions
FractionResult[4,1] <-"GR"
FractionResult[4,2] <- nrow(overlapAll)/nrow(All)
FractionResult[4,3] <- nrow(overlapMotif)/nrow(Motif)
FractionResult[4,4] <- nrow(overlapModel)/nrow(Model)
FractionResult[4,5] <- nrow(overlapAntimodel)/nrow(AntiModel)

# Make a heatmap
rownames(FractionResult) <- FractionResult[,1]
pheatmap(FractionResult[,c(2:5)], cluster_cols=F, cluster_rows=F, col=colorRampPalette(c("white","lightblue","blue","black"))(100), breaks = seq(0,1,length.out = 100))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 3](../Links/Figure3.md)