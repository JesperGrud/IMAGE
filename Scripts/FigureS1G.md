```R
# Source the necessary libaries
library(org.Hs.eg.db)
library(clusterProfiler)
library(KEGG.db)
library(goseq)

# Import the data
load("Data/3T3.R")

## GO analysis
Prediction <- Target_Genes$PPARG_1.motif
Convert <- bitr(Prediction$Factor, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
Prediction <- merge(Prediction, Convert, by.x="Factor", by.y="SYMBOL")
Prediction <- merge(Prediction, GeneRPKM[,c("countLength","Factor")], by="Factor")
KEGG <- as.list(KEGGPATHID2EXTID)
KEGG <- KEGG[grep("hsa",names(KEGG))]
Names <- as.data.frame(KEGGPATHID2NAME)
Names$path_id <- paste("hsa", Names$path_id, sep="")

Target.Vector <- as.integer(Prediction[,"ENTREZID"] %in% Prediction[ Prediction$Target == 1, "ENTREZID"])
Length.Vector <- as.integer(Prediction$countLength)
names(Target.Vector) <- Prediction$ENTREZID
names(Length.Vector) <- Prediction$ENTREZID

Target.nullp <- nullp(Target.Vector, genome="hg19", id="refGene",bias.data=Length.Vector)
Target.EA <- goseq(Target.nullp, genome="hg19",id="refGene", gene2cat=KEGG, use_genes_without_cat = T)
Target.EA <- merge(Target.EA, Names, by.x="category", by.y="path_id", sort=F)
Target.EA <- Target.EA[ Target.EA$numInCat <= 1000 & Target.EA$numInCat >= 20 & Target.EA$numDEInCat >= 1,]
Target.EA$FDR <- p.adjust(Target.EA$over_represented_pvalue, method="fdr")

# Plot the significant ones (top 4)
barplot(-log10(Target.EA[c(1:4),"FDR"]), ylim=c(0,3), las=2, ylab="-log10 FDR-corrected P-value", names=Target.EA[c(1:4),"path_name"])
```

[Back to start](../README.md)<br>
[Back to overview of Figure S1](../Links/FigureS1.md)