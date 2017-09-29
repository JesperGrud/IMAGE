```R
# Import the data
load("Data/3T3.R")

# Import the site counts and peak calls
Jun <- read.table("Data/cJun.pos", quote="\"", stringsAsFactors=FALSE)
Junb <- read.table("Data/Junb.pos", quote="\"", stringsAsFactors=FALSE)
PPAR <- read.table("/Data/PPARG.pos", quote="\"", stringsAsFactors=FALSE)
GR <- read.table("Data/GR.pos", quote="\"", stringsAsFactors=FALSE)
Counts <- read.delim("Data/Site.counts")

# Get targets and boxplot tag enrichment
par(mfcol=c(1,4))
All <- Target_Sites$PPARG_1.motif
Motif <- All[ All$NMotif > 0,]
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
AntiModel <- All[ All$NMotif > 0 & All$Contribution < 0,]
boxplot(Counts[ Counts[,1] %in% All$PeakID,12],Counts[ Counts[,1] %in% Motif$PeakID,12],Counts[ Counts[,1] %in% Model$PeakID,12],Counts[ Counts[,1] %in% AntiModel$PeakID,12], outline=F, las=1, main="PPARg")
All <- Target_Sites$JUN_1.motif
Motif <- All[ All$NMotif > 0,]
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
AntiModel <- All[ All$NMotif > 0 & All$Contribution < 0,]
boxplot(Counts[ Counts[,1] %in% All$PeakID,9],Counts[ Counts[,1] %in% Motif$PeakID,9],Counts[ Counts[,1] %in% Model$PeakID,9],Counts[ Counts[,1] %in% AntiModel$PeakID,9], outline=F, las=1, main="Jun")
All <- Target_Sites$JUNB_1.motif
Motif <- All[ All$NMotif > 0,]
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
AntiModel <- All[ All$NMotif > 0 & All$Contribution < 0,]
boxplot(Counts[ Counts[,1] %in% All$PeakID,11],Counts[ Counts[,1] %in% Motif$PeakID,11],Counts[ Counts[,1] %in% Model$PeakID,11],Counts[ Counts[,1] %in% AntiModel$PeakID,11], outline=F, las=1, main="Junb")
All <- Target_Sites$NR3C1_1.motif
Motif <- All[ All$NMotif > 0,]
Model <- All[ All$NMotif > 0 & All$Contribution > 0,]
AntiModel <- All[ All$NMotif > 0 & All$Contribution < 0,]
boxplot(Counts[ Counts[,1] %in% All$PeakID,11],Counts[ Counts[,1] %in% Motif$PeakID,11],Counts[ Counts[,1] %in% Model$PeakID,11],Counts[ Counts[,1] %in% AntiModel$PeakID,11], outline=F, las=1, main="NR3C1")
```

[Back to start](../README.md)<br>
[Back to overview of Figure 3](../Links/Figure3.md)