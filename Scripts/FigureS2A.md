```R
# Import the data
load("Data/hMSC.R")

# Plot the expression patterns
par(mfcol=c(2,3))
barplot(as.matrix(GeneRPKM[ GeneRPKM$Factor == "TSHZ1",c(74:79)]), ylim=c(0,8), main="TSHZ1", las=1, names=c("D0","4h","D1","D3","D7","D14"))
barplot(as.matrix(GeneRPKM[ GeneRPKM$Factor == "MYBL1",c(74:79)]), ylim=c(0,8), main="MYBL1", las=1, names=c("D0","4h","D1","D3","D7","D14"))
barplot(as.matrix(GeneRPKM[ GeneRPKM$Factor == "MAZ",c(74:79)]), ylim=c(0,8), main="MAZ", las=1, names=c("D0","4h","D1","D3","D7","D14"))
barplot(as.matrix(GeneRPKM[ GeneRPKM$Factor == "HSF1",c(74:79)]), ylim=c(0,8), main="HSF1", las=1, names=c("D0","4h","D1","D3","D7","D14"))
barplot(as.matrix(GeneRPKM[ GeneRPKM$Factor == "NFIL3",c(74:79)]), ylim=c(0,8), main="NFIL3", las=1, names=c("D0","4h","D1","D3","D7","D14"))
barplot(as.matrix(GeneRPKM[ GeneRPKM$Factor == "SATB1",c(74:79)]), ylim=c(0,8), main="SATB1", las=1, names=c("D0","4h","D1","D3","D7","D14"))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S2](../Links/FigureS2.md)