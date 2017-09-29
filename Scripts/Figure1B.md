```R
# Source the necessary libraries
library(pheatmap)

# Read the data
Matrix <- read.delim("Data/SREBF2_correlation.txt")
rownames(Matrix) <- Matrix[,1]

# Plot the heatmap
pheatmap(Matrix[,c(2:7)], col=colorRampPalette(c("white","darkblue"))(100))
```

[Back to start](../README.md)<br>
[Back to overview of Figure 1](../Links/Figure1.md)