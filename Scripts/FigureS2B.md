```R
# Import the data
load("Data/hMSC.R")

par(mfcol=c(2,3))
plot(c(1,2,3,4,5,6),Result[ Result$Factor == "TSHZ1",c(4:9)], type="l", ylab="Motif activity", las=1, xaxt="n", xlab="Time point", main="TSHZ1")
axis(1, labels = c("D0","4h","D1","D3","D7","D14"), at=c(1,2,3,4,5,6))
plot(c(1,2,3,4,5,6),Result[ Result$Factor == "MYBL1" & Result$Motif == "MYBL1_1.motif",c(4:9)], type="l", ylab="Motif activity", las=1, xaxt="n", xlab="Time point", main="MYBL1")
axis(1, labels = c("D0","4h","D1","D3","D7","D14"), at=c(1,2,3,4,5,6))
plot(c(1,2,3,4,5,6),Result[ Result$Factor == "MAZ" & Result$Motif == "MAZ_1.motif",c(4:9)], type="l", ylab="Motif activity", las=1, xaxt="n", xlab="Time point", main="MAZ")
axis(1, labels = c("D0","4h","D1","D3","D7","D14"), at=c(1,2,3,4,5,6))
plot(c(1,2,3,4,5,6),Result[ Result$Factor == "HSF1",c(4:9)], type="l", ylab="Motif activity", las=1, xaxt="n", xlab="Time point", main="HSF1")
axis(1, labels = c("D0","4h","D1","D3","D7","D14"), at=c(1,2,3,4,5,6))
plot(c(1,2,3,4,5,6),Result[ Result$Factor == "NFIL3",c(4:9)], type="l", ylab="Motif activity", las=1, xaxt="n", xlab="Time point", main="NFIL3")
axis(1, labels = c("D0","4h","D1","D3","D7","D14"), at=c(1,2,3,4,5,6))
plot(c(1,2,3,4,5,6),Result[ Result$Factor == "SATB1",c(4:9)], type="l", ylab="Motif activity", las=1, xaxt="n", xlab="Time point", main="SATB1")
axis(1, labels = c("D0","4h","D1","D3","D7","D14"), at=c(1,2,3,4,5,6))
```

[Back to start](../README.md)<br>
[Back to overview of Figure S2](../Links/FigureS2.md)