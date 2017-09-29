```R
# Load the data
Data <- read.delim("Data/ORO_Quant.txt", header=T)

# Process the data
Processed <- data.frame(matrix(ncol=4,nrow=8))
Processed[,1] <- c("UT","Ctrl","HSF1","MYBL1","NFIL3","MAZ","SATB1","TSHZ1")
Processed[1,2] <- mean(colMeans(Data[,c(1:2)]))
Processed[2,2] <- mean(colMeans(Data[,c(3:4)]))
Processed[3,2] <- mean(colMeans(Data[,c(5:6)]))
Processed[4,2] <- mean(colMeans(Data[,c(7:8)]))
Processed[5,2] <- mean(colMeans(Data[,c(9:10)]))
Processed[6,2] <- mean(colMeans(Data[,c(11:12)]))
Processed[7,2] <- mean(colMeans(Data[,c(13:14)]))
Processed[8,2] <- mean(colMeans(Data[,c(15:16)]))
Processed[1,3] <- sd(colMeans(Data[,c(1:2)]))/sqrt(2)
Processed[2,3] <- sd(colMeans(Data[,c(3:4)]))/sqrt(2)
Processed[3,3] <- sd(colMeans(Data[,c(5:6)]))/sqrt(2)
Processed[4,3] <- sd(colMeans(Data[,c(7:8)]))/sqrt(2)
Processed[5,3] <- sd(colMeans(Data[,c(9:10)]))/sqrt(2)
Processed[6,3] <- sd(colMeans(Data[,c(11:12)]))/sqrt(2)
Processed[7,3] <- sd(colMeans(Data[,c(13:14)]))/sqrt(2)
Processed[8,3] <- sd(colMeans(Data[,c(15:16)]))/sqrt(2)
Processed[1,4] <- 1
Processed[2,4] <- t.test(c(Data[,1], Data[,2]), c(Data[,3], Data[,4]))$p.value
Processed[3,4] <- t.test(c(Data[,1], Data[,2]), c(Data[,5], Data[,6]))$p.value
Processed[4,4] <- t.test(c(Data[,1], Data[,2]), c(Data[,7], Data[,8]))$p.value
Processed[5,4] <- t.test(c(Data[,1], Data[,2]), c(Data[,9], Data[,10]))$p.value
Processed[6,4] <- t.test(c(Data[,1], Data[,2]), c(Data[,11], Data[,12]))$p.value
Processed[7,4] <- t.test(c(Data[,1], Data[,2]), c(Data[,13], Data[,14]))$p.value
Processed[8,4] <- t.test(c(Data[,1], Data[,2]), c(Data[,15], Data[,16]))$p.value
Processed[,5] <- p.adjust(Processed[,4], method = "bonferroni")

# Plot the result
B <- barplot(Processed[,2], las=2, ylab="Oil red O stained area", ylim=c(0,5), names=Processed[,1], col=ifelse(Processed[,5] <= 0.05, ifelse(Processed[,2] <= 1, "red","green"),"grey"))
arrows(B, Processed[,2]+Processed[,3], B, Processed[,2]-Processed[,3], code=3, angle = 90, length = 0.05)
```

[Back to start](../README.md)<br>
[Back to overview of Figure 5](../Links/Figure5.md)