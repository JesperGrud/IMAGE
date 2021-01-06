#!/usr/bin/Rscript

## Catch and parse arguments from command-line
args <- commandArgs(trailingOnly=TRUE)
Arguments <- as.data.frame(args)
Directory <- as.character(Arguments[1,1])
String <- as.character(Arguments[2,1])
RegionFile <- as.character(Arguments[3,1])
GeneFile <- as.character(Arguments[4,1])
Pairing <- as.numeric(as.character(Arguments[5,1]))
ResultName <- as.character(Arguments[6,1])
Processors <- as.numeric(as.character(Arguments[7,1]))
Thres <- as.numeric(as.character(Arguments[8,1]))

Arguments[,2] <- grepl("no", Arguments[,1])
Arguments[,3] <- grepl("yes", Arguments[,1])
for (i in 1:nrow(Arguments)) {
	if (Arguments[i,2] == "TRUE" | Arguments[i,3] == "TRUE") { Index <- i }
	}
TargetsOption <- as.character(Arguments[Index,1])
Stage <- as.character(Arguments[(Index+1),1])

EnhancerDesign <- as.numeric(as.character(Arguments[c(9:(Index-1)),1]))
RNADesign <- as.numeric(as.character(Arguments[((Index+2):nrow(Arguments)),1]))

## Silence warnings
options(warn=-1)

## Load required packages silently
suppressPackageStartupMessages(library(methods))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(edgeR))

#### Stage 1 of script
if (Stage == "1") {
## Combine variables for input and output files
OutputLocation <- paste(ResultName,".R", sep="")
CountFile <- paste(Directory,"/tmp/",String,".count", sep="")
ConvertFile <- paste(Directory,"/utils/Genename_Motif.txt", sep="")

## Read data
cat("\n\t\t1. Reading data into R")
Input <- read.delim(RegionFile, header=FALSE)
Genes <- read.delim(GeneFile, header=TRUE)
Convert <- read.delim(ConvertFile, header=FALSE)
colnames(Convert) <- c("Factor","Motif","Evidence")
Counts <- fread(CountFile, verbose=FALSE, showProgress=FALSE, data.table=F)

### Process RNA-seq data
## Define variables for gene expression data dimensions
Samples <- length(RNADesign)
Conditions <- length(unique(RNADesign))

## Calculate RPKM and differential expression
cat("\n\t\t2. Analyzing gene expression")
rownames(Genes) <- Genes$Symbol
GeneDGE <- DGEList(counts = Genes[,c(9:(8+Samples))], group = RNADesign)
GeneDGE <- calcNormFactors(GeneDGE)
GeneRPKM <- rpkm(GeneDGE, gene.length=Genes$countLength, log=T, normalized.lib.sizes=T)
GeneRPKM <- as.data.frame(rpkm(GeneDGE, gene.length=Genes$countLength, log=T, normalized.lib.sizes=T))
for (samp in 1:Samples) { colnames(GeneRPKM)[samp] <- paste("Expression-Sample_",samp,"-Condition_", RNADesign[samp], sep="") }
GeneRPKM$Max <- apply(GeneRPKM[,c(1:Samples)],1,FUN="max")
GeneRPKM$Factor <- rownames(GeneRPKM)
GeneRPKM <- GeneRPKM[ GeneRPKM$Max >= Thres,]

if (Conditions == 2) {
GeneDGE <- estimateCommonDisp(GeneDGE)
GeneDGE <- estimateTrendedDisp(GeneDGE)
GeneDGE <- estimateTagwiseDisp(GeneDGE)
Test <- as.data.frame(topTags(exactTest(GeneDGE),nrow(Genes)))
Test$Factor <- rownames(Test)
colnames(Test)[1] <- "logFC_Cond2-vs-Cond1"
colnames(Test)[3] <- "Pval_Cond2-vs-Cond1"
colnames(Test)[4] <- "FDR_Cond2-vs-Cond1"
GeneRPKM <- merge(GeneRPKM, Test[,c(1,3,4,5)], by="Factor")
}

if (Conditions >= 3) {
if (Pairing == 0) {
DesignMatrix <- model.matrix(~factor(RNADesign))
GeneDGE <- estimateGLMCommonDisp(GeneDGE, DesignMatrix)
GeneDGE <- estimateGLMTrendedDisp(GeneDGE, DesignMatrix)
GeneDGE <- estimateGLMTagwiseDisp(GeneDGE, DesignMatrix)
GeneFit <- glmQLFit(GeneDGE, DesignMatrix)

Tests <- expand.grid(RNADesign, RNADesign)
Tests <- Tests[ Tests[,1] != Tests[,2],]
Tests$Combine <- paste(Tests[,1], Tests[,2], sep="-")
Tests <- Tests[duplicated(Tests$Combine)==F,c(1,2)]
Tests <- Tests[ Tests[,1] < Tests[,2],]
Tests <- Tests[ order(Tests[,1], Tests[,2]),]

Coefficients <- Tests[ Tests[,1] == 1,2]
Contrasts <- Tests[ Tests[,1] > 1,]

for (i in 1:nrow(Contrasts)) {
Contrast <- as.vector(rep(0,Conditions))
Contrast[Contrasts[i,1]] <- -1
Contrast[Contrasts[i,2]] <- 1
Test <- as.data.frame(topTags(glmQLFTest(GeneFit, contrast=Contrast),nrow(Genes)))
Test$Factor <- rownames(Test)
colnames(Test)[1] <- paste("logFC_Cond",Contrasts[i,2],"-vs-Cond",Contrasts[i,1], sep="")
colnames(Test)[4] <- paste("Pval_Cond",Contrasts[i,2],"-vs-Cond",Contrasts[i,1], sep="")
colnames(Test)[5] <- paste("FDR_Cond",Contrasts[i,2],"-vs-Cond",Contrasts[i,1], sep="")
GeneRPKM <- merge(GeneRPKM, Test[,c(1,4,5,6)], by="Factor")
}

for (i in Coefficients) {
Test <- as.data.frame(topTags(glmQLFTest(GeneFit, coef=i),nrow(Genes)))
Test$Factor <- rownames(Test)
colnames(Test)[1] <- paste("logFC_Cond",i,"-vs-Cond1", sep="")
colnames(Test)[4] <- paste("Pval_Cond",i,"-vs-Cond1", sep="")
colnames(Test)[5] <- paste("FDR_Cond",i,"-vs-Cond1", sep="")
GeneRPKM <- merge(GeneRPKM, Test[,c(1,4,5,6)], by="Factor")
}
}

if (Pairing == 1) {
Table <- table(RNADesign)
Paired <- vector()
for (q in 1:length(Table)) { Paired <- c(Paired, seq(1,Table[q], by=1))}

DesignMatrix <- model.matrix(~factor(Paired) + factor(RNADesign))
GeneDGE <- estimateGLMCommonDisp(GeneDGE, DesignMatrix)
GeneDGE <- estimateGLMTrendedDisp(GeneDGE, DesignMatrix)
GeneDGE <- estimateGLMTagwiseDisp(GeneDGE, DesignMatrix)
GeneFit <- glmQLFit(GeneDGE, DesignMatrix)

Tests <- expand.grid(RNADesign, RNADesign)
Tests <- Tests[ Tests[,1] != Tests[,2],]
Tests$Combine <- paste(Tests[,1], Tests[,2], sep="-")
Tests <- Tests[duplicated(Tests$Combine)==F,c(1,2)]
Tests <- Tests[ Tests[,1] < Tests[,2],]
Tests <- Tests[ order(Tests[,1], Tests[,2]),]

MaxNumberOfReplicates <- max(table(RNADesign))

Coefficients <- Tests[ Tests[,1] == 1,2] + (MaxNumberOfReplicates - 1)
Contrasts <- Tests[ Tests[,1] > 1,]
Contrasts[,3] <- Contrasts[,1] + (MaxNumberOfReplicates - 1)
Contrasts[,4] <- Contrasts[,2] + (MaxNumberOfReplicates - 1)

for (i in 1:nrow(Contrasts)) {
Contrast <- as.vector(rep(0,(Conditions+(MaxNumberOfReplicates - 1))))
Contrast[Contrasts[i,3]] <- -1
Contrast[Contrasts[i,4]] <- 1
Test <- as.data.frame(topTags(glmQLFTest(GeneFit, contrast=Contrast),nrow(Genes)))
Test$Factor <- rownames(Test)
colnames(Test)[1] <- paste("logFC_Cond",Contrasts[i,2],"-vs-Cond",Contrasts[i,1], sep="")
colnames(Test)[4] <- paste("Pval_Cond",Contrasts[i,2],"-vs-Cond",Contrasts[i,1], sep="")
colnames(Test)[5] <- paste("FDR_Cond",Contrasts[i,2],"-vs-Cond",Contrasts[i,1], sep="")
GeneRPKM <- merge(GeneRPKM, Test[,c(1,4,5,6)], by="Factor")
}

for (i in Coefficients) {
Test <- as.data.frame(topTags(glmQLFTest(GeneFit, coef=i),nrow(Genes)))
Test$Factor <- rownames(Test)
colnames(Test)[1] <- paste("logFC_Cond",i,"-vs-Cond1", sep="")
colnames(Test)[4] <- paste("Pval_Cond",i,"-vs-Cond1", sep="")
colnames(Test)[5] <- paste("FDR_Cond",i,"-vs-Cond1", sep="")
GeneRPKM <- merge(GeneRPKM, Test[,c(1,4,5,6)], by="Factor")
}
}
}

# Merge results and set gene names to capitals
GeneRPKM <- merge(GeneRPKM, Genes[,c(1:8)], by.x="Factor", by.y="Symbol")
GeneRPKM$Factor <- toupper(GeneRPKM$Factor )

### Process motif data
## Set the column names on the matrix of motif occurences
cat("\n\t\t3. Converting motif hits to matrix - Takes a while")
colnames(Counts)[1] <- "ID"
colnames(Counts)[4] <- "Motif"

## Remove motifs that pair with TFs expressed below threshold
Expressed <- Convert
Expressed$Expression <- 0
Expressed[ Expressed$Factor %in% GeneRPKM$Factor,"Expression"] <- 1
Expressed <- Expressed[ order(Expressed$Motif, -Expressed$Expression),]
Expressed <- Expressed[ duplicated(Expressed$Motif) == F,]
Expressed <- Expressed[ Expressed$Expression == 1,]
Counts <- Counts[ Counts$Motif %in% Expressed$Motif,]

## Convert list of motif hits into count matrix
Counts$Count <- 1
Motifs <- Counts[ duplicated(Counts$Motif)==F,"Motif"]
CountMatrix <- data.frame(ID = Input$V4)
for (i in 1:length(Motifs)) {
	CurrentMotif <- Motifs[i]
	Tmp <- Counts[ Counts[,4] == CurrentMotif,]
	Tmp2 <- aggregate(Tmp$Count, by=list(Tmp$ID), FUN="sum")
	colnames(Tmp2) <- c("ID",as.character(Tmp[1,4]))
	CountMatrix <- merge(CountMatrix, Tmp2, all.x=T, by="ID")
	CountMatrix[ is.na(CountMatrix[,(i+1)]), (i+1)] <- 0
}

## Order by PeakID, save the PeakIDs as rownames and remove the column
rownames(CountMatrix) <- CountMatrix[,1]
CountMatrix <- CountMatrix[ order(CountMatrix[,1]),]
CountMatrix <- CountMatrix[,c(2:ncol(CountMatrix))]

## Define variable with the number of motifs
Motifs <- ncol(CountMatrix)

## Keep a copy of CountMatrix for target site and target gene predictions
CountMatrixUnnorm <- CountMatrix
CountMatrixUnnorm$ID <- rownames(CountMatrixUnnorm)
## Normalize the count matrix
List <- data.frame(ID=rownames(CountMatrix))
CountMatrix <- as.matrix(CountMatrix)
Cols <- colMeans(CountMatrix)
CountMatrix <- t(t(CountMatrix) - Cols)
CountMatrix <- as.data.frame(CountMatrix)
CountMatrix <- cbind(List, CountMatrix)

### Process the ChIP-seq data
## Define variables for ChIP-seq data dimensions
Samples <- length(EnhancerDesign)
Conditions <- length(unique(EnhancerDesign))

## Set column and rownames for chip-seq data
colnames(Input)[4] <- "ID" 
rownames(Input) <- Input[,4]

## Normalize the matrices
cat("\n\t\t4. Calculating motif activity - Full model stage")
Input <- Input[ order(Input[,4]), c(5:ncol(Input))]
Input <- as.matrix(Input)
Cols <- colMeans(Input)
SD <- apply(Input[,c(1:ncol(Input))],2,FUN="sd")
Input <- t( ((t(Input) - Cols) / SD) ) 
Input <- as.data.frame(Input)
Input$ID <- rownames(Input)
Input <- Input[,c(ncol(Input),1:(ncol(Input)-1))]

### Perform ridge regression of all samples
cat("\n\t\t5. Performing ridge regression")
## Setup a list to catch all regressions
RegressionList <- list()
## Register the needed amount of workers 
registerDoParallel(Processors)
## Define the objects needed to catch to output
Coef <- matrix(ncol=Samples, nrow=Motifs)

## Loop through samples
for (samp in 1:Samples) {
	# Setup data.frame for testing
	Test <- merge(Input[,c(1,(1+samp))], CountMatrix, by="ID")
	Test <- Test[,c(2:ncol(Test))]
	colnames(Test)[1] <- "Exprs"
	# Determine lambda using 10-fold cross validation in parallel (glmnet package)
	cv <- cv.glmnet(y=as.matrix(Test[,1]), x=as.matrix(Test[,c(2:ncol(Test))]), alpha=0, standardize=FALSE, intercept=TRUE, parallel=TRUE, standardize.response=FALSE)
	lambda <- cv$lambda.min
	# Do the regression
	out <- glmnet(y=as.matrix(Test[,1]), x=as.matrix(Test[,c(2:ncol(Test))]), alpha=0, standardize=FALSE, intercept=TRUE, lambda=lambda, standardize.response=FALSE)
	# Capture entire regression in a list
	RegressionList[[samp]] <- out
	# Capture the regression coefficients
	Tmp <- as.data.frame(as.matrix(coef(out)))
	Coef[,samp] <- Tmp[c(2:nrow(Tmp)),1]
	# Print the progress
	cat(sprintf("\n\t\t\tSample %s completed", samp))
	Sys.sleep(5)
	flush.console()
}

### Process the ridge regression data
## Setup objects with regression data 
cat("\n\t\t6. Calculating motif activity - Reduced model stage")
Coef <- data.frame(Coef)
colnames(Coef) <- EnhancerDesign
Coef$Motif <- colnames(CountMatrix)[c(2:ncol(CountMatrix))]
Sys.sleep(5)
flush.console()
## Save a copy of the normalized CountMatrix for later use
Backup <- CountMatrix
## Rearrange CountMatrix to make summation easy
CountMatrix <- t(CountMatrix)
CountMatrix <- data.frame(CountMatrix)
CountMatrix$Motif <- rownames(CountMatrix)
CountMatrix <- CountMatrix[c(2:nrow(CountMatrix)),]
CountMatrix <- CountMatrix[ order(CountMatrix$Motif),]
## Calculate Eps - Npm*Ams and summarize pr sample
FullModel <- matrix(nrow=(ncol(CountMatrix)-1), ncol=Samples)
for (samp in 1:Samples) {
	Tmp2 <- matrix(ncol=(ncol(CountMatrix)-1), nrow=nrow(CountMatrix))
	Tmp <- Coef[,c(samp,(Samples+1))]
	Tmp <- Tmp[ order(Tmp$Motif),]
	for (i in 1:(ncol(CountMatrix)-1)) { Tmp2[,i] <- as.numeric(as.character(CountMatrix[,i])) * Tmp[,1] }
	FullModel[,samp] <- (as.numeric(as.character(Input[,(1+samp)])) - colSums(Tmp2))^2
}
## Calculate average standard deviation
Average <- sum(rowSums(FullModel))/(nrow(FullModel)*ncol(FullModel))
## Split R session to reduce memory impact of parallelization
save(list = c("Coef","CountMatrix","Samples","Input","FullModel","Average"), file = paste(Directory,"/tmp/",String,"_Parallel.tmp", sep=""), envir = .GlobalEnv)
save(list = setdiff(ls(),c("Coef","CountMatrix","Samples","Input","FullModel","Average")), file = paste(Directory,"/tmp/",String,"_TheRest.tmp", sep=""), envir = .GlobalEnv)
}

#### Stage 2 of script
if (Stage == "2") {
## Load only what is needed for parallel processing
load(paste(Directory,"/tmp/",String,"_Parallel.tmp", sep=""))
## Initialize needed spaces
RedModel <- matrix(nrow=(ncol(CountMatrix)-1), ncol=Samples)
Tmp2 <- matrix(ncol=(ncol(CountMatrix)-1), nrow=nrow(CountMatrix)-1)
Score <- matrix(nrow=nrow(RedModel), ncol=1)
CoefList <- list()
## Setup to do parallel processing
registerDoParallel(cores = Processors)
## Perform target prediction in parallel (!HIGH MEMORY USAGE)
Targets <- foreach(motf=1:nrow(Coef), .inorder=TRUE) %dopar% {
	## Set stringsAsFactors to false to prevent empty output
	options(stringsAsFactors = FALSE)
	## Define the motif being investigated
	Motif <- Coef[motf,"Motif"]
	## Remove the motif from the count matrix and make a new count matrix
	TestCountMatrix <- CountMatrix[ rownames(CountMatrix) != Motif,]
	TestCountMatrix <- TestCountMatrix[ order(TestCountMatrix$Motif),]
	## Convert new matrix to matrix structure
	for (i in 1:(ncol(TestCountMatrix)-1)) { Tmp2[,i] <- as.numeric(as.character(TestCountMatrix[,i])) }
	## Get coefficients and place into a list for all samples
	for (samp in 1:Samples) {
		Tmp <- Coef[,c(samp,(Samples+1))]
		Tmp <- Tmp[ Tmp$Motif != Motif,]
		Tmp <- Tmp[ order(Tmp$Motif),]
		CoefList[[samp]] <- Tmp
	}
	## Calculate Eps - Npm*Ams and summarize pr sample
	for (samp in 1:Samples) {
		Out <- Tmp2 * as.numeric(as.character(CoefList[[samp]][,1]))
		RedModel[,samp] <- (as.numeric(Input[,(1+samp)]) - colSums(Out))^2
	}
	## Calculate the score
	for (site in 1:nrow(RedModel)) { Score[site,1] <- sum(RedModel[site,] - FullModel[site,])/Average }
	## Put results into a target list
	return(Score)
}

## Load the rest of the working space
load(paste(Directory,"/tmp/",String,"_TheRest.tmp", sep=""))
### Process the predicted target sites and enhancer activities
Regions <- read.delim(RegionFile, header=FALSE)
colnames(Regions)[4] <- "ID"
Target_Sites <- Targets
rm(Targets)

## Set names and statistics of site predictions
for (motf in 1:nrow(Coef)) {
	## Define the motif being investigated
	Motif <- Coef[motf,"Motif"]
	## Extract the PeakIDs and motif counts
	Matrix <- CountMatrixUnnorm[ , colnames(CountMatrixUnnorm) == Motif | colnames(CountMatrixUnnorm) == "ID"]
	## Extract the scores
	Tmp <- Target_Sites[[motf]]
	## Combine the results
	Tmp <- cbind(Matrix, Tmp)
	## Merge with enhancer information, rearrange and set column names
	Tmp <- merge(Tmp, Regions[,c(1:4)], by="ID")
	Tmp <- Tmp[,c(4,5,6,1,2,3)]
	colnames(Tmp) <- c("Chr","Start","End","PeakID","NMotif","Contribution")
	## Call target sites
	Tmp$Target <- 0
	Tmp[ Tmp[,5] > 0 & Tmp[,6] > 0, "Target"] <- 1
	## Reinsert the new data.frame into the list
	Target_Sites[[motf]] <- Tmp
	names(Target_Sites)[[motf]] <- Motif
}

## Process enhancer activity
Enhancer_Activity <- Coef
for (samp in 1:Samples) { colnames(Enhancer_Activity)[samp] <- paste("Activity-Sample_",samp,"-Condition_", EnhancerDesign[samp], sep="") }

### Process the predicted enhancers to setup for stage 2 regression analysis
cat("\n\t\t7. Creating weight matrices")
## Define TSS for all genes
Tmp <- GeneRPKM[,c("Chr","Start","End","Strand","Factor")]
colnames(Tmp)[c(2:3)] <- c("TSS","TSS")
TSS <- rbind(Tmp[ Tmp[,4] == "-", c(1,3,5)], Tmp[ Tmp[,4] == "+", c(1,2,5)])
TSS$Start <- TSS$TSS - 100000
TSS$End <- TSS$TSS + 100000
TSS$subjectHits <- seq(1, nrow(TSS), by=1)
GenesRanges <- GRanges(seqnames=TSS[,1], IRanges(start=TSS$Start, end=TSS$End), strand=rep("+", nrow(TSS)), names=TSS[,3])

## Make a GRange object from the enhancers
Regions$queryHits <- seq(1, nrow(Regions), by=1)
RegionRanges <- GRanges(seqnames=Regions[,1], IRanges(start=Regions[,2], end=Regions[,3]), strand=rep("+", nrow(Regions)), names=Regions[,4])

## Overlap gene regions and enhancers
Overlap <- findOverlaps(query=RegionRanges, subject=GenesRanges, type="any")
Overlap <- as.data.frame(Overlap)

## Merge the appropriate information from gene and enhancers
Overlap <- merge(Overlap, Regions[,c(1:4,ncol(Regions))], by="queryHits")
Overlap <- merge(Overlap, TSS[,c("Factor","TSS","subjectHits")], by="subjectHits")

## Calculate scaled regulatory potential
Overlap$Center <- (Overlap[,5]+Overlap[,4])/2
Overlap$Distance <- Overlap$Center - Overlap$TSS
Overlap <- Overlap[ abs(Overlap$Distance) <= 100000,]
Overlap$Potential <- (exp(-(0.5+4*(abs(Overlap$Distance)/100000)))-exp(-1 * (0.5 + 4)))/(max((exp(-1 * (0.5 + (4 * (seq(0,100000,by=10)/100000))))-exp(-1 * (0.5 + 4)))))

## Extract the motif count at predicted target sites
## OBS! The WeightMatrix is not reproducible from the original IMAGE pipeline. Problem must be localized in this loop, which has been redefined.
Scores <- data.frame(matrix(nrow=nrow(Target_Sites[[1]]), ncol=nrow(Coef)))
for (q in 1:nrow(Coef)) { 
	Tmp2 <- Target_Sites[[q]]
	Tmp2[ Tmp2[,6] < 0,5] <- 0
	colnames(Tmp2)[5] <- names(Target_Sites)[q]
	Tmp2 <- Tmp2[ order(Tmp2[,4]),]
	Scores[,q] <- Tmp2[,5]
}
## Combine the results
colnames(Scores) <- names(Target_Sites)
Tmp2 <- Target_Sites[[q]]
Tmp2 <- Tmp2[ order(Tmp2[,4]),]
Scores <- cbind(Tmp2$PeakID, Scores)
colnames(Scores)[1] <- "ID"
Overlap <- merge(Overlap, Scores, by="ID")

## Multiply the motif count with the distance weight
for (q in 1:nrow(Coef)) { Overlap[,(q+11)] <- Overlap[,(q+11)] * Overlap[,11] }
	
## Summarize for each gene
WeightMatrix <- aggregate(Overlap[,c(12:(nrow(Coef)+11))], by=list(Overlap$Factor), FUN="sum")
colnames(WeightMatrix)[1] <- "Factor"

## Normalize the motif matrix
List <- data.frame(Factor=WeightMatrix$Factor)
WeightMatrix <- WeightMatrix[,c(2:ncol(WeightMatrix))]
WeightMatrix <- as.matrix(WeightMatrix)
Cols <- colMeans(WeightMatrix)
WeightMatrix <- t(t(WeightMatrix) - Cols)
WeightMatrix <- as.data.frame(WeightMatrix)
WeightMatrix <- cbind(List, WeightMatrix)

### Process RNA-seq data to prepare for stage 2 regression
## Define variables for gene expression data dimensions
Samples <- length(RNADesign)
Conditions <- length(unique(RNADesign))

## Normalize the gene expression matrix by Z-transformation
rownames(GeneRPKM) <- GeneRPKM$Factor
Input <- GeneRPKM[,c(2:(Samples+1))]
Input <- as.matrix(Input)
Cols <- colMeans(Input)
SD <- apply(Input[,c(1:ncol(Input))],2,FUN="sd")
Input <- t( ((t(Input) - Cols) / SD) ) 
Input <- as.data.frame(Input)
Input$Factor <- rownames(Input)
Input <- Input[,c(ncol(Input),1:(ncol(Input)-1))]

### Perform stage 2 ridge regression of all samples
cat("\n\t\t8. Calculating factor activity - Full model stage")
## Setup a list to catch all regressions
RegressionList <- list()
## Register the needed amount of workers 
registerDoParallel(Processors)
## Define the objects needed to catch to output
Coef <- matrix(ncol=Samples, nrow=Motifs)
## Loop through samples
for (samp in 1:Samples) {
	# Setup data.frame for testing
	Test <- merge(Input[,c(1,(1+samp))], WeightMatrix, by="Factor")
	Test <- Test[,c(2:ncol(Test))]
	colnames(Test)[1] <- "Exprs"
	# Determine lambda using 10-fold cross validation in parallel (glmnet package)
	cv <- cv.glmnet(y=as.matrix(Test[,1]), x=as.matrix(Test[,c(2:ncol(Test))]), alpha=0, standardize=FALSE, intercept=TRUE, parallel=TRUE, standardize.response=FALSE)
	lambda <- cv$lambda.min
	# Do the regression
	out <- glmnet(y=as.matrix(Test[,1]), x=as.matrix(Test[,c(2:ncol(Test))]), alpha=0, standardize=FALSE, intercept=TRUE, lambda=lambda, standardize.response=FALSE)
	# Capture entire regression in a list
	RegressionList[[samp]] <- out
	# Capture the regression coefficients
	Tmp <- as.data.frame(as.matrix(coef(out)))
	Coef[,samp] <- Tmp[c(2:nrow(Tmp)),1]
	# Print the progress
	cat(sprintf("\n\t\t\tSample %s completed", samp))
	Sys.sleep(5)
	flush.console()
}

## Setup objects
Coef <- data.frame(Coef)
colnames(Coef) <- RNADesign
Coef$Motif <- colnames(WeightMatrix)[2:(Motifs+1)]

## Convert to standard scores for significance testing
Significant <- Coef
for (i in 1:Samples) {
SD <- sd(Significant[,i])
Mean <- mean(Significant[,i])
Significant[,i] <- (Significant[,i] - Mean) / SD
}

## Set up for all combinations of tests
Tests <- expand.grid(RNADesign, RNADesign)
Tests <- Tests[ Tests[,1] != Tests[,2],]
Tests$Combine <- paste(Tests[,1], Tests[,2], sep="-")
Tests <- Tests[duplicated(Tests$Combine)==F,c(1,2)]
Tests <- Tests[ Tests[,1] < Tests[,2],]
Tests <- Tests[ order(Tests[,1], Tests[,2]),]	
SampleOutline <- cbind(RNADesign, seq(1, Samples, by=1))

## Perform T tests
for (q in 1:nrow(Tests)) {
Significant[,(ncol(Significant)+1)] <- 1
colnames(Significant)[ncol(Significant)] <- paste("Pvalue_", Tests[q,1] ,"_vs_", Tests[q,2],  sep="")
if (Pairing == 1) {
for (i in 1:nrow(Significant)) { Significant[i,ncol(Significant)] <- t.test(as.numeric(Significant[i,SampleOutline[ SampleOutline[,1] == Tests[ q,1],2]]),as.numeric(Significant[i,SampleOutline[ SampleOutline[,1] == Tests[ q,2],2]]), paired=TRUE)$p.value }
} else {
for (i in 1:nrow(Significant)) { Significant[i,ncol(Significant)] <- t.test(as.numeric(Significant[i,SampleOutline[ SampleOutline[,1] == Tests[ q,1],2]]),as.numeric(Significant[i,SampleOutline[ SampleOutline[,1] == Tests[ q,2],2]]))$p.value }
}
}


## Get the minimal P-values for each motif
if (length(grep("Pvalue", colnames(Significant))) == 1) {
Significant$Min <- Significant[,grep("Pvalue", colnames(Significant))]
}

if (length(grep("Pvalue", colnames(Significant))) >= 2) {
Significant$Min <- apply(Significant[,grep("Pvalue", colnames(Significant))],1,FUN="min")
}

## Identify putative regulators by motif activity and gene expression
if (length(grep("FDR", colnames(GeneRPKM))) == 1) {
GeneRPKM$MinimalFDR <- GeneRPKM[,grep("FDR", colnames(GeneRPKM))]
}

if (length(grep("FDR", colnames(GeneRPKM))) >= 2) {
GeneRPKM$MinimalFDR <- apply(GeneRPKM[,grep("FDR", colnames(GeneRPKM))],1,FUN="min")
}

Significant <- merge(Significant, Convert, by="Motif")
Significant <- merge(Significant, GeneRPKM[,c("MinimalFDR","Factor","Max")], all.x=T, by="Factor")
Significant[ is.na(Significant$MinimalFDR), "MinimalFDR"] <- 1
Significant[ is.na(Significant$Max), "Max"] <- 1

## Calculate average motif activity for each sample
ColIndices <- data.frame(Design = as.numeric(RNADesign), Col = seq(3,(2+length(RNADesign)),by=1))
for (i in unique(RNADesign)) {
Significant$Tmp <- 1
Significant[,ncol(Significant)] <- rowMeans(Significant[,ColIndices[ ColIndices$Design == i, 2]])
colnames(Significant)[ncol(Significant)] <- paste("Activity_Condition_", i, sep="")
}

## Get maximal activity
Significant$MaximalActivity <- apply(Significant[ ,grep("Activity_Condition", colnames(Significant))],1,FUN="max")

## Calculate the weighted P-value
Significant$Combined <- (rank(-Significant$Max) / nrow(Significant)) * (rank(-Significant$MaximalActivity) / nrow(Significant)) * (rank(Significant$Min) / nrow(Significant)) * (rank(Significant$MinimalFDR) / nrow(Significant))
Full <- Significant
Significant <- Full[ Full$Combined <= 0.001,]
Significant <- Significant[ Significant$Min <= 0.05,]
Significant <- Significant[ duplicated(Significant$Motif) ==F,]

### Process target genes
cat("\n\t\t9. Calculating factor activity - Reduced model stage")
Sys.sleep(5)
flush.console()
## Calculate Eps - Wpm*Ams and summarize pr sample
Input <- Input[ Input$Factor %in% WeightMatrix$Factor,]
Input <- Input[ order(Input$Factor),]
Tmp3 <- WeightMatrix
Tmp3 <- Tmp3[ order(Tmp3[,1]),]
Tmp3 <- Tmp3[,c(2:ncol(Tmp3))]
Tmp3 <- t(as.matrix(Tmp3))
Tmp3 <- as.data.frame(Tmp3)
Tmp3$Motif <- rownames(Tmp3)
Tmp3 <- Tmp3[ order(Tmp3$Motif),]
FullModel <- matrix(nrow=(ncol(Tmp3)-1), ncol=Samples)
Tmp4 <- Coef
Tmp4 <- Tmp4[ duplicated(Tmp4$Motif) == F,]
for (samp in 1:Samples) {
	Tmp2 <- matrix(ncol=(ncol(Tmp3)-1), nrow=nrow(Tmp3))
	Tmp <- Tmp4[,c((samp),ncol(Tmp4))]
	Tmp <- Tmp[ order(Tmp$Motif),]
	for (i in 1:(ncol(Tmp3)-1)) { Tmp2[,i] <- as.numeric(as.character(Tmp3[,i])) * Tmp[,1] }
	FullModel[,samp] <- (as.numeric(as.character(Input[,(1+samp)])) - colSums(Tmp2))^2
}
## Calculate average standard deviation
Average <- sum(rowSums(FullModel))/(nrow(FullModel)*ncol(FullModel))
## Define the motifs to test - dependent on the target variable
TestMotifs <- Coef
if (TargetsOption == "no") {
	TestMotifs <- TestMotifs[ TestMotifs$Motif %in% Significant$Motif,]
}
## Split R session to reduce memory impact of parallelization
save(list = c("TestMotifs","WeightMatrix","Samples","Input","FullModel","Average", "Tmp4"), file = paste(Directory,"/tmp/",String,"_Parallel_2.tmp", sep=""), envir = .GlobalEnv)
save(list = setdiff(ls(),c("TestMotifs","WeightMatrix","Samples","Input","FullModel","Average","Tmp4")), file = paste(Directory,"/tmp/",String,"_TheRest_2.tmp", sep=""), envir = .GlobalEnv)
}

#### Stage 3 of script
if (Stage == "3") {
## Load only what is needed for parallel processing
load(paste(Directory,"/tmp/",String,"_Parallel_2.tmp", sep=""))
## Setup to do parallel processing
registerDoParallel(Processors)
## Loop through all motifs and determine target genes
Targets <- foreach(motf=1:nrow(TestMotifs), .inorder=TRUE) %dopar% {
	## Set stringsAsFactors to false to prevent empty output
	options(stringsAsFactors = FALSE)
	## Define the motif being investigated
	Motif <- TestMotifs[motf,"Motif"]
	## Define a matrix to capture the reduced model
	RedModel <- matrix(nrow=(nrow(WeightMatrix)), ncol=Samples)
	## Remove the motif from the weight matrices
	Tmp3 <- WeightMatrix
	Tmp3 <- Tmp3[ ,colnames(Tmp3) != Motif]
	Tmp3 <- Tmp3[ order(Tmp3[,1]),]
	Tmp3 <- Tmp3[,c(2:ncol(Tmp3))]
	Tmp3 <- t(as.matrix(Tmp3))
	Tmp3 <- as.data.frame(Tmp3)
	Tmp3$Motif <- rownames(Tmp3)
	Tmp3 <- Tmp3[ order(Tmp3$Motif),]
	for (samp in 1:Samples) {
		Tmp2 <- matrix(ncol=(ncol(Tmp3)-1), nrow=nrow(Tmp3))
		Tmp <- Tmp4[,c((samp),ncol(Tmp4))]
		Tmp <- Tmp[ Tmp$Motif != Motif,]
		Tmp <- Tmp[ order(Tmp$Motif),]
		for (i in 1:(ncol(Tmp3)-1)) { Tmp2[,i] <- as.numeric(as.character(Tmp3[,i])) * Tmp[,1] }
		RedModel[,samp] <- (as.numeric(as.character(Input[,(1+samp)])) - colSums(Tmp2))^2
	}
	## Define a matrix for catching the score
	Score <- matrix(nrow=nrow(RedModel), ncol=1)
	## Calculate the score
	for (site in 1:nrow(RedModel)) { Score[site,1] <- sum(RedModel[site,]- FullModel[site,])/Average }
	## Calculate the maximal summary of weights
	Matrix <- data.frame(WeightMatrix)
	Matrix <- Matrix[, colnames(Matrix) == Motif | colnames(Matrix) == "Factor"]
	## Combine the score and occurences
	Score <- cbind(Matrix, Score)
	## Calculate a p-value based on rank
	Score$Pval <- (rank(-Score[,2]) / nrow(Score)) * (rank(-Score[,3]) / nrow(Score))
	## Put results into a target list
	return(Score)
}

## Load the reset of the working space
load(paste(Directory,"/tmp/",String,"_TheRest_2.tmp", sep=""))

### Process the results and finalize output
## Define putative hits
cat("\n\t\t10. Processing results")

## Setting up all activities and enhancer regions for saving
Gene_Activity <- Coef
Enhancers <- Regions[,c(1:(ncol(Regions)-1))]

## Processing target genes, setting names and statistics
for (motf in 1:nrow(TestMotifs)) {
	Motif <- TestMotifs[motf,"Motif"]
	Tmp <- Targets[[motf]]
	Tmp$Target <- 0
	Tmp[ Tmp$Pval <= 0.005, "Target"] <- 1
	Targets[[motf]] <- Tmp
	names(Targets)[[motf]] <- Motif
}
Target_Genes <- Targets
## Calculate average motif activity for each sample, and remove replicate data
ColIndices <- data.frame(Design = as.numeric(RNADesign), Col = seq(1,(length(RNADesign)),by=1))
for (i in unique(RNADesign)) {
Coef$Tmp <- 1
Coef[,ncol(Coef)] <- rowMeans(Coef[,ColIndices[ ColIndices$Design == i, 2]])
colnames(Coef)[ncol(Coef)] <- paste("Activity_Condition_", i, sep="")
}
## Merge information
Coef <- suppressWarnings(merge(Coef, Full[,c(1,2,grep("Evidence",colnames(Full)),grep("Min",colnames(Full)),grep("Combined",colnames(Full)))], by="Motif"))
Coef <- Coef[ , grep("MinimalFDR", colnames(Coef), invert=T)]
## Calculate average gene expression strength
ColIndices <- data.frame(Design = as.numeric(RNADesign), Col = seq(2,(1+length(RNADesign)),by=1))
for (i in unique(RNADesign)) {
GeneRPKM$Tmp <- 1
GeneRPKM[,ncol(GeneRPKM)] <- rowMeans(GeneRPKM[,ColIndices[ ColIndices$Design == i, 2]])
colnames(GeneRPKM)[ncol(GeneRPKM)] <- paste("Expression_Condition_", i, sep="")
}
## Perform correlations between motif activity and gene expression (per condition)
Coef <- merge(Coef, GeneRPKM[,c(1,grep("Expression-Sample",colnames(GeneRPKM)))], by="Factor")
Coef$Pearsons <- 0
for (i in 1:nrow(Coef)) { Coef[i,"Pearsons"] <- cor(as.numeric(Coef[i,c(3:(2+Samples))]),as.numeric(Coef[i,grep("Expression-Sample",colnames(Coef))])) }

## Setup and set column names
Coef <- Coef[,grep("Expression-Sample", colnames(Coef), invert=T)]
Coef <- Coef[,c(1,2,(Samples+3):ncol(Coef))]
colnames(Coef)[grep("Min",colnames(Coef))] <- "MinimalActivityPvalue"
colnames(Coef)[grep("Combined",colnames(Coef))] <- "WeightedPvalue"

## Identify putatitative regulatory factors
Coef$CausalTF <- 0
Coef[ Coef$WeightedPvalue <= 0.005 & Coef$MinimalActivityPvalue <= 0.05, "CausalTF"] <- 2
Coef[ Coef$WeightedPvalue <= 0.001 & Coef$MinimalActivityPvalue <= 0.05, "CausalTF"] <- 1

Coef <- Coef[,c(1,2,grep("Evidence",colnames(Coef)),3:(2+Conditions),(grep("Evidence",colnames(Coef))+1):ncol(Coef))]
Result <- Coef
rm(list=setdiff(ls(), c("Target_Sites","Target_Genes","GeneRPKM","OutputLocation","Enhancer_Activity","Gene_Activity","Enhancers","Result")))
cat("\n\t\t11. Outputting results\n\n")
# Save an R session image
save.image(OutputLocation)
}
