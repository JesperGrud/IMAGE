```R
# Load the necessary libraries
library(pheatmap)

# Import and merge data
ReadIn <- read.csv("Data/HSF_metabolites.txt", sep="\t", header=T)
colnames(ReadIn)[c(2)] <- "HSF_Flux"
Data <- ReadIn[,c(3,4,5,6,14,17,2)]

ReadIn <- read.csv("Data/MYBL1.metabolites.txt", sep="\t", header=T)
colnames(ReadIn)[c(2)] <- "MYBL1_Flux"
Data <- merge(Data, ReadIn[,c(2,3)], by="Reaction.Abbreviation")

ReadIn <- read.csv("Data/MAZ.metabolites.txt", sep="\t", header=T)
colnames(ReadIn)[c(2)] <- "MAZ_Flux"
Data <- merge(Data, ReadIn[,c(2,3)], by="Reaction.Abbreviation")

ReadIn <- read.csv("Data/SATB1.metabolites.txt", sep="\t", header=T)
colnames(ReadIn)[c(2)] <- "SATB1_Flux"
Data <- merge(Data, ReadIn[,c(2,3)], by="Reaction.Abbreviation")

ReadIn <- read.csv("Data/TSHZ1.metabolites.txt", sep="\t", header=T)
colnames(ReadIn)[c(2)] <- "TSHZ1_Flux"
Data <- merge(Data, ReadIn[,c(2,3)], by="Reaction.Abbreviation")

ReadIn <- read.csv("Data/Uni_D0.metabolites.txt", sep="\t", header=T)
colnames(ReadIn)[c(2)] <- "Uni_D0_Flux"
Data <- merge(Data, ReadIn[,c(2,3)], by="Reaction.Abbreviation")

ReadIn <- read.csv("Data/Uni_D7.metabolites.txt", sep="\t", header=T)
colnames(ReadIn)[c(2)] <- "Uni_D7_Flux"
Data <- merge(Data, ReadIn[,c(2,3)], by="Reaction.Abbreviation")
Data <- Data[ Data$Reaction.Name != "",]

# Calculate maximal flux
Data$Max <- apply(Data[,c(7:11,13)],1,FUN="max")

# Find changed metabolites with flux above 1e-4 and plot a heatmap
Data <- Data[ !is.na(Data$HSF_Flux),]
rownames(Data) <- Data$Reaction.Abbreviation
pheatmap(Data[ (Data$Subsystem == "Triacylglycerol synthesis" | Data$Subsystem == "Fatty acid synthesis") & Data$Max >= 1e-04 & Data$Confidence.Level > 0,c(12,13,7:11)], cluster_cols = F, scale="row")
```

[Back to start](../README.md)<br>
[Back to overview of Figure S2](../Links/FigureS2.md)