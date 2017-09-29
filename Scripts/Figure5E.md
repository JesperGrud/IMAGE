```R
# Import libraries
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(goseq)
library(reactome.db)

# Set the working directory
setwd("C:/Data/Projects/PhD/Bioinformatics/Analyser/Motifs/Final/Analysis")

# Clean the working space
rm(list=ls())

# Import the data
Counts <- read.delim("Data/Candidates.count.txt")
rownames(Counts) <- Counts$RefSeq

# Setup to analyze using edgeR
dd <- data.frame(a = gl(2,2))
Design <- model.matrix(~ a, dd)

# Analyze using edgeR
DGE <- DGEList(Counts[,c(45,46,49,50)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
Control <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
Control$RefSeq <- rownames(Control)

DGE <- DGEList(Counts[,c(49,50,13,14)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
HSF <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
HSF$RefSeq <- rownames(HSF)

DGE <- DGEList(Counts[,c(49,50,19,20)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
MAZ <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
MAZ$RefSeq <- rownames(MAZ)

DGE <- DGEList(Counts[,c(49,50,25,26)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
MYBL1 <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
MYBL1$RefSeq <- rownames(MYBL1)

DGE <- DGEList(Counts[,c(49,50,37,38)])
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
SATB1 <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
SATB1$RefSeq <- rownames(SATB1)

DGE <- DGEList(Counts[,c(49,50,43,44)], group = c("TSHZ1","TSHZ1","Control","Control"))
DGE <- calcNormFactors(DGE)
DGE <- estimateGLMCommonDisp(DGE, design = Design)
DGE <- estimateGLMTrendedDisp(DGE, design = Design)
DGE <- estimateGLMTagwiseDisp(DGE, design = Design)
Fit <- glmQLFit(DGE, design = Design)
TSHZ1 <- as.data.frame(topTags(glmQLFTest(Fit, coef = 2), n = nrow(Counts)))
TSHZ1$RefSeq <- rownames(TSHZ1)

## Pathway analysis
# Setup the reactome database
Reactome <- as.data.frame(reactomeEXTID2PATHID)
Relation <- read.delim("Data/ReactomePathwaysRelation.txt", header=FALSE)
Pathways <- read.delim("Data/ReactomePathways.txt", header=FALSE)
Pathways <- Pathways[ Pathways$V3 == "Homo sapiens",]
Pathways$DB_ID <- substr(Pathways$V1, 7, nchar(as.character(Pathways$V1)))

Metabolism <- Relation[ Relation[,1] %in% as.character(Pathways[ Pathways$V2 == "Metabolism",1]),]  
Metabolism <- rbind(Metabolism, Relation[ Relation[,1] %in% Metabolism[,2],])
Metabolism <- Metabolism[ duplicated(Metabolism[,2])==F,]

Current <- nrow(Metabolism)
New <- Current + 1
while (New > Current) {
  Current <- nrow(Metabolism)
  Metabolism <- rbind(Metabolism, Relation[ Relation[,1] %in% Metabolism[,2],])
  Metabolism <- Metabolism[ duplicated(Metabolism[,2])==F,]
  New <- nrow(Metabolism)
}

Pathways <- Pathways[ Pathways$V1 %in% Metabolism[,1] | Pathways$V1 %in% Metabolism[,2],]
Pathways <- Pathways[ duplicated(Pathways$V1)==F,]
Reactome <- Reactome[ Reactome$DB_ID %in% Pathways$V1,]
Reactome <- split(Reactome$gene_id, f = Reactome$DB_ID, drop=T)

# Convert to ENTREZ IDs
Convert <- bitr(Counts[ , "RefSeq"], "REFSEQ", "ENTREZID", OrgDb = org.Hs.eg.db)
Convert <- Convert[ duplicated(Convert$REFSEQ)==F,]
Convert <- Convert[ duplicated(Convert$ENTREZID)==F,]
Control <- merge(Control, Convert, by.x="RefSeq", by.y="REFSEQ")
Control <- merge(Control, Counts[,c("RefSeq","countLength")], by="RefSeq")
HSF <- merge(HSF, Convert, by.x="RefSeq", by.y="REFSEQ")
MAZ <- merge(MAZ, Convert, by.x="RefSeq", by.y="REFSEQ")
TSHZ1 <- merge(TSHZ1, Convert, by.x="RefSeq", by.y="REFSEQ")
MYBL1 <- merge(MYBL1, Convert, by.x="RefSeq", by.y="REFSEQ")
SATB1 <- merge(SATB1, Convert, by.x="RefSeq", by.y="REFSEQ")

# Make vector objects
Control_Up.Vector <- as.integer(Control[,"ENTREZID"] %in% Control[ Control[,2] > 0 & Control[,5] <= 0.05, "ENTREZID"])
Control_Down.Vector <- as.integer(Control[,"ENTREZID"] %in% Control[ Control[,2] < 0 & Control[,5] <= 0.05, "ENTREZID"])
HSF_Less.Vector <- as.integer(HSF[,"ENTREZID"] %in% HSF[ HSF[,2] <= 0 & HSF[,5] <= 0.05, "ENTREZID"])
MAZ_Less.Vector <- as.integer(MAZ[,"ENTREZID"] %in% MAZ[ MAZ[,2] <= 0 & MAZ[,5] <= 0.05, "ENTREZID"])
TSHZ1_Less.Vector <- as.integer(TSHZ1[,"ENTREZID"] %in% TSHZ1[ TSHZ1[,2] <= 0 & TSHZ1[,5] <= 0.05, "ENTREZID"])
MYBL1_Less.Vector <- as.integer(MYBL1[,"ENTREZID"] %in% MYBL1[ MYBL1[,2] <= 0 & MYBL1[,5] <= 0.05, "ENTREZID"])
SATB1_Less.Vector <- as.integer(SATB1[,"ENTREZID"] %in% SATB1[ SATB1[,2] <= 0 & SATB1[,5] <= 0.05, "ENTREZID"])
HSF_More.Vector <- as.integer(HSF[,"ENTREZID"] %in% HSF[ HSF[,2] >= 0 & HSF[,5] <= 0.05, "ENTREZID"])
MAZ_More.Vector <- as.integer(MAZ[,"ENTREZID"] %in% MAZ[ MAZ[,2] >= 0 & MAZ[,5] <= 0.05, "ENTREZID"])
TSHZ1_More.Vector <- as.integer(TSHZ1[,"ENTREZID"] %in% TSHZ1[ TSHZ1[,2] >= 0 & TSHZ1[,5] <= 0.05, "ENTREZID"])
MYBL1_More.Vector <- as.integer(MYBL1[,"ENTREZID"] %in% MYBL1[ MYBL1[,2] >= 0 & MYBL1[,5] <= 0.05, "ENTREZID"])
SATB1_More.Vector <- as.integer(SATB1[,"ENTREZID"] %in% SATB1[ SATB1[,2] >= 0 & SATB1[,5] <= 0.05, "ENTREZID"])
Length.Vector <- as.integer(Control[, "countLength"])
names(Control_Up.Vector) <- Control[,"ENTREZID"]
names(Control_Down.Vector) <- Control[,"ENTREZID"]
names(HSF_Less.Vector) <- HSF[,"ENTREZID"]
names(MAZ_Less.Vector) <- MAZ[,"ENTREZID"]
names(TSHZ1_Less.Vector) <- TSHZ1[,"ENTREZID"]
names(MYBL1_Less.Vector) <- MYBL1[,"ENTREZID"]
names(SATB1_Less.Vector) <- SATB1[,"ENTREZID"]
names(HSF_More.Vector) <- HSF[,"ENTREZID"]
names(MAZ_More.Vector) <- MAZ[,"ENTREZID"]
names(TSHZ1_More.Vector) <- TSHZ1[,"ENTREZID"]
names(MYBL1_More.Vector) <- MYBL1[,"ENTREZID"]
names(SATB1_More.Vector) <- SATB1[,"ENTREZID"]
names(Length.Vector) <- Control[,"ENTREZID"]

# Run the nullp estimation
Control_Up.nullp <- nullp(Control_Up.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
Control_Down.nullp <- nullp(Control_Down.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
HSF_Less.nullp <- nullp(HSF_Less.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
MAZ_Less.nullp <- nullp(MAZ_Less.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
TSHZ1_Less.nullp <- nullp(TSHZ1_Less.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
MYBL1_Less.nullp <- nullp(MYBL1_Less.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
SATB1_Less.nullp <- nullp(SATB1_Less.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
HSF_More.nullp <- nullp(HSF_More.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
MAZ_More.nullp <- nullp(MAZ_More.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
TSHZ1_More.nullp <- nullp(TSHZ1_More.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
MYBL1_More.nullp <- nullp(MYBL1_More.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
SATB1_More.nullp <- nullp(SATB1_More.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)

# Run enrichment analysis
Control_Up.EA <- goseq(Control_Up.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
Control_Down.EA <- goseq(Control_Down.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
HSF_Less.EA <- goseq(HSF_Less.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
MAZ_Less.EA <- goseq(MAZ_Less.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
TSHZ1_Less.EA <- goseq(TSHZ1_Less.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
MYBL1_Less.EA <- goseq(MYBL1_Less.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
SATB1_Less.EA <- goseq(SATB1_Less.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
HSF_More.EA <- goseq(HSF_More.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
MAZ_More.EA <- goseq(MAZ_More.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
TSHZ1_More.EA <- goseq(TSHZ1_More.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
MYBL1_More.EA <- goseq(MYBL1_More.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)
SATB1_More.EA <- goseq(SATB1_More.nullp, genome="hg19",id="refGene", gene2cat=Reactome, use_genes_without_cat = T)

# Clean up the results
Control_Up.EA <- Control_Up.EA[ Control_Up.EA$numInCat <= 1000 & Control_Up.EA$numInCat >= 20 & Control_Up.EA$numDEInCat >= 1,]
Control_Down.EA <- Control_Down.EA[ Control_Down.EA$numInCat <= 1000 & Control_Down.EA$numInCat >= 20 & Control_Down.EA$numDEInCat >= 1,]
HSF_Less.EA <- HSF_Less.EA[ HSF_Less.EA$numInCat <= 1000 & HSF_Less.EA$numInCat >= 20 & HSF_Less.EA$numDEInCat >= 1,]
MAZ_Less.EA <- MAZ_Less.EA[ MAZ_Less.EA$numInCat <= 1000 & MAZ_Less.EA$numInCat >= 20 & MAZ_Less.EA$numDEInCat >= 1,]
TSHZ1_Less.EA <- TSHZ1_Less.EA[ TSHZ1_Less.EA$numInCat <= 1000 & TSHZ1_Less.EA$numInCat >= 20 & TSHZ1_Less.EA$numDEInCat >= 1,]
MYBL1_Less.EA <- MYBL1_Less.EA[ MYBL1_Less.EA$numInCat <= 1000 & MYBL1_Less.EA$numInCat >= 20 & MYBL1_Less.EA$numDEInCat >= 1,]
SATB1_Less.EA <- SATB1_Less.EA[ SATB1_Less.EA$numInCat <= 1000 & SATB1_Less.EA$numInCat >= 20 & SATB1_Less.EA$numDEInCat >= 1,]
HSF_More.EA <- HSF_More.EA[ HSF_More.EA$numInCat <= 1000 & HSF_More.EA$numInCat >= 20 & HSF_More.EA$numDEInCat >= 1,]
MAZ_More.EA <- MAZ_More.EA[ MAZ_More.EA$numInCat <= 1000 & MAZ_More.EA$numInCat >= 20 & MAZ_More.EA$numDEInCat >= 1,]
TSHZ1_More.EA <- TSHZ1_More.EA[ TSHZ1_More.EA$numInCat <= 1000 & TSHZ1_More.EA$numInCat >= 20 & TSHZ1_More.EA$numDEInCat >= 1,]
MYBL1_More.EA <- MYBL1_More.EA[ MYBL1_More.EA$numInCat <= 1000 & MYBL1_More.EA$numInCat >= 20 & MYBL1_More.EA$numDEInCat >= 1,]
SATB1_More.EA <- SATB1_More.EA[ SATB1_More.EA$numInCat <= 1000 & SATB1_More.EA$numInCat >= 20 & SATB1_More.EA$numDEInCat >= 1,]

# Perform FDR-corrections
Control_Up.EA$Control_Up <- p.adjust(Control_Up.EA$over_represented_pvalue, method="fdr")
Control_Down.EA$Control_Down <- p.adjust(Control_Down.EA$over_represented_pvalue, method="fdr")
HSF_Less.EA$HSF_Less <- p.adjust(HSF_Less.EA$over_represented_pvalue, method="fdr")
MAZ_Less.EA$MAZ_Less <- p.adjust(MAZ_Less.EA$over_represented_pvalue, method="fdr")
TSHZ1_Less.EA$TSHZ1_Less <- p.adjust(TSHZ1_Less.EA$over_represented_pvalue, method="fdr")
MYBL1_Less.EA$MYBL1_Less <- p.adjust(MYBL1_Less.EA$over_represented_pvalue, method="fdr")
SATB1_Less.EA$SATB1_Less <- p.adjust(SATB1_Less.EA$over_represented_pvalue, method="fdr")
HSF_More.EA$HSF_More <- p.adjust(HSF_More.EA$over_represented_pvalue, method="fdr")
MAZ_More.EA$MAZ_More <- p.adjust(MAZ_More.EA$over_represented_pvalue, method="fdr")
TSHZ1_More.EA$TSHZ1_More <- p.adjust(TSHZ1_More.EA$over_represented_pvalue, method="fdr")
MYBL1_More.EA$MYBL1_More <- p.adjust(MYBL1_More.EA$over_represented_pvalue, method="fdr")
SATB1_More.EA$SATB1_More <- p.adjust(SATB1_More.EA$over_represented_pvalue, method="fdr")

# Merge the results
Enrichment <- merge( Pathways[,c(1,2)], Control_Up.EA[,c(1,6)], by.y="category", by.x="V1", all=T)
colnames(Enrichment)[c(1,2)] <- c("category","name") 
Enrichment <- merge(Enrichment, Control_Down.EA[,c(1,6)], by="category", all=T)                    
Enrichment <- merge(Enrichment, HSF_Less.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, HSF_More.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, MAZ_Less.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, MAZ_More.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, TSHZ1_Less.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, TSHZ1_More.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, MYBL1_Less.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, MYBL1_More.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, SATB1_Less.EA[,c(1,6)], by="category", all=T)
Enrichment <- merge(Enrichment, SATB1_More.EA[,c(1,6)], by="category", all=T)

# Set Q-value NAs to 1 
for (i in c(3:14)) { Enrichment[ is.na(Enrichment[,i]),i] <- 1 }

# Calculate enrichment in each pathway
Result <- data.frame(matrix(ncol=11, nrow=254))
for (q in 1:254) {
  Result[q,1] <- names(Reactome[q])
  A <- nrow(MYBL1[ MYBL1$FDR <= 0.05 & MYBL1$logFC > 0 & MYBL1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  B <- nrow(MYBL1[ MYBL1$FDR <= 0.05 & MYBL1$logFC < 0 & MYBL1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  C <- nrow(MYBL1[ MYBL1$FDR <= 0.05 & MYBL1$logFC > 0,])
  D <- nrow(MYBL1[ MYBL1$FDR <= 0.05 & MYBL1$logFC < 0,])
  E <- nrow(MYBL1[ MYBL1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  F <- nrow(MYBL1)
  # Add 1/1000 of the group size as pseudocount to avoid 0's
  A <- A + (C/1000)
  C <- C + (C/1000)
  B <- B + (D/1000)
  D <- D + (D/1000)
  Result[q,2] <- (A/E)/(C/F)
  Result[q,3] <- (B/E)/(D/F)
  A <- nrow(HSF[ HSF$FDR <= 0.05 & HSF$logFC > 0 & HSF$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  B <- nrow(HSF[ HSF$FDR <= 0.05 & HSF$logFC < 0 & HSF$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  C <- nrow(HSF[ HSF$FDR <= 0.05 & HSF$logFC > 0,])
  D <- nrow(HSF[ HSF$FDR <= 0.05 & HSF$logFC < 0,])
  E <- nrow(HSF[ HSF$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  F <- nrow(HSF)
  # Add 1/1000 of the group size as pseudocount to avoid 0's
  A <- A + (C/1000)
  C <- C + (C/1000)
  B <- B + (D/1000)
  D <- D + (D/1000)
  Result[q,4] <- (A/E)/(C/F)
  Result[q,5] <- (B/E)/(D/F)
  A <- nrow(SATB1[ SATB1$FDR <= 0.05 & SATB1$logFC > 0 & SATB1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  B <- nrow(SATB1[ SATB1$FDR <= 0.05 & SATB1$logFC < 0 & SATB1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  C <- nrow(SATB1[ SATB1$FDR <= 0.05 & SATB1$logFC > 0,])
  D <- nrow(SATB1[ SATB1$FDR <= 0.05 & SATB1$logFC < 0,])
  E <- nrow(SATB1[ SATB1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  F <- nrow(SATB1)
  # Add 1/1000 of the group size as pseudocount to avoid 0's
  A <- A + (C/1000)
  C <- C + (C/1000)
  B <- B + (D/1000)
  D <- D + (D/1000)
  Result[q,6] <- (A/E)/(C/F)
  Result[q,7] <- (B/E)/(D/F)
  A <- nrow(TSHZ1[ TSHZ1$FDR <= 0.05 & TSHZ1$logFC > 0 & TSHZ1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  B <- nrow(TSHZ1[ TSHZ1$FDR <= 0.05 & TSHZ1$logFC < 0 & TSHZ1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  C <- nrow(TSHZ1[ TSHZ1$FDR <= 0.05 & TSHZ1$logFC > 0,])
  D <- nrow(TSHZ1[ TSHZ1$FDR <= 0.05 & TSHZ1$logFC < 0,])
  E <- nrow(TSHZ1[ TSHZ1$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  F <- nrow(TSHZ1)
  # Add 1/1000 of the group size as pseudocount to avoid 0's
  A <- A + (C/1000)
  C <- C + (C/1000)
  B <- B + (D/1000)
  D <- D + (D/1000)
  Result[q,8] <- (A/E)/(C/F)
  Result[q,9] <- (B/E)/(D/F)
  A <- nrow(MAZ[ MAZ$FDR <= 0.05 & MAZ$logFC > 0 & MAZ$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  B <- nrow(MAZ[ MAZ$FDR <= 0.05 & MAZ$logFC < 0 & MAZ$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  C <- nrow(MAZ[ MAZ$FDR <= 0.05 & MAZ$logFC > 0,])
  D <- nrow(MAZ[ MAZ$FDR <= 0.05 & MAZ$logFC < 0,])
  E <- nrow(MAZ[ MAZ$ENTREZID %in% as.data.frame(Reactome[q])[,1],])
  F <- nrow(MAZ)
  # Add 1/1000 of the group size as pseudocount to avoid 0's
  A <- A + (C/1000)
  C <- C + (C/1000)
  B <- B + (D/1000)
  D <- D + (D/1000)
  Result[q,10] <- (A/E)/(C/F)
  Result[q,11] <- (B/E)/(D/F)
}

colnames(Result) <- c("category","MYBL1_Higher_Enrichment","MYBL1_Lower_Enrichment","HSF1_Higher_Enrichment","HSF1_Lower_Enrichment","SATB1_Higher_Enrichment","SATB1_Lower_Enrichment","TSHZ1_Higher_Enrichment","TSHZ1_Lower_Enrichment","MAZ_Higher_Enrichment","MAZ_Lower_Enrichment")
Enrichment <- merge(Enrichment[,c(1:14)], Result, by="category")

# Setup correct order
Enrichment <- Enrichment[,c(1:14,17,18,15,16,23,24,19,20,21,22)]

# Plot the pathways of interest
par(mfcol=c(1,1))
barplot(as.matrix(log2(Enrichment[ Enrichment$name == "Metabolism of lipids",c(15:24)])), las=2, main="Metabolism of lipids", ylab="log2 Enrichment")
```

[Back to start](../README.md)<br>
[Back to overview of Figure 5](../Links/Figure5.md)