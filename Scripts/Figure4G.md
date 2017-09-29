```R
# Import the data
load("Data/3T3.R")

# Test D0 vs D7
Gene_Activity$D7_Pvalue <- 1
Gene_Activity$D0_Pvalue <- 1
Gene_Activity$Pvalue <- 1
for (i in 1:nrow(Gene_Activity)) { Gene_Activity[i,"D7_Pvalue"] <- t.test(as.numeric(Gene_Activity[i,c(7,8)]))$p.value}
for (i in 1:nrow(Gene_Activity)) { Gene_Activity[i,"D0_Pvalue"] <- t.test(as.numeric(Gene_Activity[i,c(1,2)]))$p.value}
for (i in 1:nrow(Gene_Activity)) { Gene_Activity[i,"Pvalue"] <- t.test(as.numeric(Gene_Activity[i,c(1,2)]),as.numeric(Gene_Activity[i,c(7,8)]))$p.value}
Gene_Activity$D7 <- rowMeans(Gene_Activity[,c(7,8)])
Gene_Activity$D0 <- rowMeans(Gene_Activity[,c(1,2)])

# Import activators and repressors
Activators <- read.delim("Data/Activators")
Activators$Gene.names <- substr(Activators$Entry.name, 0, regexpr("_HUMAN",Activators$Entry.name)-1)
Repressors <- read.delim("Data/Repressors")
Repressors$Gene.names <- substr(Repressors$Entry.name, 0, regexpr("_HUMAN",Repressors$Entry.name)-1)

Repressed <- Gene_Activity[ Gene_Activity$D7 < -0.005 & Gene_Activity$D7_Pvalue <= 0.05 ,]
Induced <- Gene_Activity[ Gene_Activity$D7 > 0.005 & Gene_Activity$D7_Pvalue <= 0.05,]

# Import all TFs
TFs <- read.delim("Data/Genename_Motif.txt", header=FALSE)

# Calculate jaccard index (Correlation treshold = 0.7) and take the log2 ratio
X <- nrow(TFs[ TFs[,1] %in% Repressors$Gene.names,])
Y <- nrow(TFs[ TFs[,1] %in% Activators$Gene.names,])

Tmp <- TFs[ TFs[,2] %in% Repressed[ ,"Motif"] & TFs[,1] %in% Repressors$Gene.names & !(TFs[,1] %in% Activators$Gene.names),]
Tmp <- Tmp[ duplicated(Tmp$V1)==F,]
A <- nrow(Tmp)
Tmp <- TFs[ TFs[,2] %in% Repressed[ ,"Motif"] & TFs[,1] %in% Activators$Gene.names & !(TFs[,1] %in% Repressors$Gene.names),]
Tmp <- Tmp[ duplicated(Tmp$V1)==F,]
B <- nrow(Tmp)
Tmp <- TFs[ TFs[,2] %in% Repressed[ ,"Motif"],]
Tmp <- Tmp[ duplicated(Tmp$V1)==F,]
C <- nrow(Tmp)
Repressed <- log2((B / (C + Y -  B)) / (A / (C + X -  A)))

Tmp <- TFs[ TFs[,2] %in% Induced[ ,"Motif"] & TFs[,1] %in% Repressors$Gene.names & !(TFs[,1] %in% Activators$Gene.names),]
Tmp <- Tmp[ duplicated(Tmp$V1)==F,]
A <- nrow(Tmp)
Tmp <- TFs[ TFs[,2] %in% Induced[ ,"Motif"] & TFs[,1] %in% Activators$Gene.names & !(TFs[,1] %in% Repressors$Gene.names),]
Tmp <- Tmp[ duplicated(Tmp$V1)==F,]
B <- nrow(Tmp)
Tmp <- TFs[ TFs[,2] %in% Induced[ ,"Motif"],]
Tmp <- Tmp[ duplicated(Tmp$V1)==F,]
C <- nrow(Tmp)
Induced <- log2((B / (C + Y -  B)) / (A / (C + X -  A)))

# Plot the results
barplot(c(Induced, Repressed), las=1, ylab="Activator to repressor log2 ratio", col=c("green","red"), names=c("Pos","Neg"), ylim=c(-1,1))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)