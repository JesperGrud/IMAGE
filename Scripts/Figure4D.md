```R
# Source the necessary libraries
library(pheatmap)

# Import the data
load("Data/3T3.R")

# Calculate the mean expression of the predicted target genes
Expression <- data.frame(matrix(ncol=5, nrow=726))
for (i in 1:726) {
  Targets <- Target_Genes[[i]]
  Targets <- Targets[ Targets$Target == 1,]
  Expression[i,1] <- names(Target_Genes)[i]
  Expression[i,c(2:5)] <- colMeans(t(scale(t(GeneRPKM[ GeneRPKM$Factor %in% Targets$Factor,c(37,38,39,40)]))))
}

# Merge the results
Result <- merge(Result, Expression, by.x="Motif", by.y="X1")
Result <- Result[ Result$CausalTF == 1,]

# Plot a heatmap of the scaled motif activities
Result$ID <- paste(Result$Factor, Result$Motif, sep="-")
rownames(Result) <- Result$ID
P <- pheatmap(Result[ ,c(4:7)], cluster_rows = T, cluster_cols = F, border_color = NA, scale="row", clustering_distance_rows = "maximum", color = colorRampPalette(c("red","white","green"))(100))
pheatmap(Result[ ,c(12,13,14,15)], cluster_rows = P$tree_row, cluster_cols = F, border_color = NA, scale="row", color = colorRampPalette(c("red","white","green"))(100))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 4](../Links/Figure4.md)