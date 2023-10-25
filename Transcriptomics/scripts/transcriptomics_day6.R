# Class notes and script from October 25, 2023
# Carolyn Hanrahan 

## Set your working directory
setwd("C:/Users/Carolyn/Desktop/Ecological Genomics/Ecological_Genomics/Transcriptomics/data")

# Load the package

if (!require("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")

BiocManager::install("GO.db")


library(WGCNA);

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

library(DESeq2)
library(ggplot2)

library(tidyverse)

library(remotes)
remotes::install_github("kevinblighe/CorLevelPlot")

library(CorLevelPlot) 
library(gridExtra)

library(Rmisc) 

#A. Import the counts matrix and metadata and filter using DESeq2
# 1. Import the counts matrix and metadata and filter using DESeq2

countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)