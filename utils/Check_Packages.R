#!/usr/bin/Rscript

if (is.element('glmnet', installed.packages()[,1]) == F) {
cat("FALSE")
}

if (is.element('methods', installed.packages()[,1]) == F) {
cat("FALSE")
}

if (is.element('cluster', installed.packages()[,1]) == F) {
cat("FALSE")
}

if (is.element('doParallel', installed.packages()[,1]) == F) {
cat("FALSE")
}

if (is.element('foreach', installed.packages()[,1]) == F) {
cat("FALSE")
}

if (is.element('Matrix', installed.packages()[,1]) == F) {
cat("FALSE")
}

if (is.element('data.table', installed.packages()[,1]) == F) {
cat("FALSE")
}

if (is.element('GenomicRanges', installed.packages()[,1]) == F) {
cat("FALSE")
}
