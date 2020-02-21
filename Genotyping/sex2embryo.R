args <- commandArgs(TRUE)

# libraries
library(cellrangerRkit)
library(GenomicFeatures)
require(ggplot2)

set.seed(42)

# data
cellranger_dir <- args[1]
cell2embryo_output <- args[2]
outDir <- args[3]
ID <- args[4]

# files
out_RData <- file.path(outDir, paste(ID, "_sex2embryo.RData", sep=""))

# variables
genome <- "mm10"

XY <- c("ENSMUSG00000096768", "ENSMUSG00000069045", "ENSMUSG00000069049")
XX <- "ENSMUSG00000086503"


# load cellranger output and extract RNA expression information
cellranger <- load_cellranger_matrix(cellranger_dir, genome=genome)
expr <- as.matrix(exprs(cellranger))

# reduce to XX and XY determining genes
expr <- expr[c(XX,XY), ]

# load BC (cell) to cluster (embryo) matrix
load(cell2embryo_output)

if(!(length(grep("-1", colnames(expr)))>0 & length(grep("-1", names(cluster_assignment)))>0)){
    names(cluster_assignment) <- paste0(names(cluster_assignment), "-1")
}

# calculate fraction of cells per cluster expressing XX or XY genes
sex_classification <- data.frame(cluster=integer(), XX=integer(), XY=integer())
for(i in 1:max(cluster_assignment)){
    XXexpr <- mean(expr[XX, which(colnames(expr) %in% names(cluster_assignment[cluster_assignment==i]))]>0)
    XYexpr <- mean(expr[XY, which(colnames(expr) %in% names(cluster_assignment[cluster_assignment==i]))]>0)
    sex_classification[i,] <- c(i, XXexpr, XYexpr)
}
save(centers, cluster_assignment, sex_classification, XX, XY, file=out_RData)
