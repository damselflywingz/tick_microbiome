## ---------------------------
##
## Script name: 2022_11_07_heatmap.R
##
## Purpose of script: Generates a heatmap of log2CPM expression from meta-transcriptome data
##
## Author: Dr. Amber R Paulson
##
## Date Created: 2022-11-07
##
## 
## Email: Amber.Rose.Paulson@gmail.com
## Github: https://github.com/damselflywingz
##
## ---------------------------
##
## Notes:
## See the README.md for more information on how to get started - https://github.com/damselflywingz/tick_microbiome
## !IMPORTANT - The R Project should be opened with RStudio first, and then next proceed to open and run the script.
##   
##
## ---------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

## read in a table containg the raw count data for the meta-transcriptome 
geneExp <- read.table("rawGeneExpression.txt", header=T,row.names = 1, sep="\t")
geneExp <- geneExp[1:4]

BiocManager::install("limma")
BiocManager::install("statmod")
library("limma")
library("statmod")

## filter non-expressed genes (at least raw count of 5 in 1 or more libraries)

filter <- apply(geneExp, 1, function(x) length(x[x>4])>=1)
myfiltered <- geneExp[filter, ]

BiocManager::install("edgeR")
library(edgeR)

myxx <- as.factor(c("38F","38M","42F","45F"))

## create the DGEList
countsDGEList <- DGEList(counts=myfiltered, group=myxx)

## normalize
UQdge <- calcNormFactors(countsDGEList, method="upperquartile")

library(limma)

design <- model.matrix(~0+myxx)
myvm <- voom(UQdge, design=design, plot=TRUE)

## see https://rstudio-pubs-static.s3.amazonaws.com/378564_de40904f6df2462cbf8ae7134e99d9e8.html

install.packages("tidyr")
library(tidyr)

## organize the log2CPM expression data for the heatmap
vm.norm <- myvm$E
vm.norm2 <- cbind(rownames(vm.norm),data.frame(vm.norm,row.names=NULL))
colnames(vm.norm2) <- c("GeneName","X38F","X38M","X42F","X45F")
vm.norm2$GeneName <- factor(vm.norm2$GeneName, levels = vm.norm2$GeneName)
df_sp = vm.norm2 %>% gather(sample, value, -GeneName)

library(ggplot2)
library(viridis)

ggplot(df_sp, aes(x = sample, y = GeneName, fill = value)) + geom_tile() + scale_fill_viridis(option="plasma") + ylab("GeneName") + xlab("samples")
