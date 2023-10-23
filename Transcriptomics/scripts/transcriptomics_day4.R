# Class Notes: ahud_DESeq2
# Carolyn Hanrahan 
# October 18 2023 


## Set your working directory
setwd("C:/Users/Carolyn/Desktop/Ecological Genomics/Ecological_Genomics/Transcriptomics/data")

## Import the libraries that we're likely to need in this session

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

BiocManager::install("vsn")
BiocManager::install("hexbin")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ read in data ~~~~~~~~~~~~


# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample discription table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ explore data distribution ~~~~~~


# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,21000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # [1] 8217.81
median(rowSums(countsTableRound)) # [1] 377

apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,1000), ylim=c(0,50000),breaks=10000)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ working with DESeq, filtering ~~~~

#### Create a DESeq object and define the experimental design here with the tilda

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ generation + treatment)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 30) >= 28,]
nrow(dds) 

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

###############################################################

# Check the quality of the data by sample clustering and visualization
# The goal of transformation "is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low."

library("pheatmap")
library("vsn")

# this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))


sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
###############################################################

# PCA to visualize global gene expression patterns

# first transform the data for plotting using variance stabilization
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



############################################################### MORE ADVANCED PCA plots

# Let's plot the PCA by generation in four panels

data <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

###########  

dataF0 <- subset(data, generation == 'F0')

F0 <- ggplot(dataF0, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-10, 25) + xlim(-40, 10)+ # zoom for F0 with new assembly
  #ylim(-40, 25) + xlim(-50, 50)+ # new assembly limits
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  ##theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

F0


#png("PCA_F0.png", res=300, height=5, width=5, units="in")

#ggarrange(F0, nrow = 1, ncol=1)

#dev.off()

################# F2

dataF2 <- subset(data, generation == 'F2')

F2 <- ggplot(dataF2, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23), labels = c("Ambient", "Acidification","Warming"))+
  # scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) + 
  #scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A"), labels = c("Ambient", "Acidification","Warming"))+
  theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #scale_size(guide="none") +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F2


# png("PCA_F2.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F2, nrow = 1, ncol=1)
# 
# dev.off()

# Yes - F2 is missing one ambient replicate

################################ F4

dataF4 <- subset(data, generation == 'F4')

F4 <- ggplot(dataF4, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+ # limits with filtered assembly
  #ylim(-20, 10) + xlim(-40, 25)+  # zoom with filtered assembly
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  # scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) + 
  #scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  #theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #scale_size(guide="none") +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F4


# png("PCA_F4.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F4, nrow = 1, ncol=1)
# 
# dev.off()


################# F11

dataF11 <- subset(data, generation == 'F11')

F11 <- ggplot(dataF11, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-45, 25) + xlim(-50, 55)+
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,24), labels = c("Ambient", "OWA"))+
  scale_fill_manual(values=c('#6699CC', "#CC3333"), labels = c("Ambient", "OWA"))+
  guides(shape = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F11


# png("PCA_F11.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F11, nrow = 1, ncol=1)
# 
# dev.off()

ggarrange(F0, F2, F4, F11, nrow = 2, ncol=2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Differential Expresion

## Check on the DE results from the DESeq command way above from resultsNames(dds)

resAM_OWA <- results(dds, name="treatment_OWA_vs_AM", alpha=0.05)

resAM_OWA <- resAM_OWA[order(resAM_OWA$padj),]
head(resAM_OWA)  

summary(resAM_OWA)


resAM_OW <- results(dds, name="treatment_OW_vs_AM", alpha=0.05)

resAM_OW <- resAM_OW[order(resAM_OW$padj),] #ordered by adjusted p-value, displays out most significant reads 

# Results 
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747             470.7031
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.744::m.744            3876.8668
# TRINITY_DN1275_c0_g1::TRINITY_DN1275_c0_g1_i9::g.6875::m.6875       169.5286
# TRINITY_DN17372_c0_g1::TRINITY_DN17372_c0_g1_i11::g.43161::m.43161   87.1611
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733            698.5519
# TRINITY_DN1272_c0_g1::TRINITY_DN1272_c0_g1_i7::g.6644::m.6644       572.5323
# log2FoldChange
# <numeric>
#   TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747                  1.750149
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.744::m.744                  1.424818
# TRINITY_DN1275_c0_g1::TRINITY_DN1275_c0_g1_i9::g.6875::m.6875           -1.646946
# TRINITY_DN17372_c0_g1::TRINITY_DN17372_c0_g1_i11::g.43161::m.43161       2.503328
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733                 2.203755
# TRINITY_DN1272_c0_g1::TRINITY_DN1272_c0_g1_i7::g.6644::m.6644           -0.782739
# lfcSE
# <numeric>
#   TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747             0.305583
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.744::m.744             0.260806
# TRINITY_DN1275_c0_g1::TRINITY_DN1275_c0_g1_i9::g.6875::m.6875       0.343458
# TRINITY_DN17372_c0_g1::TRINITY_DN17372_c0_g1_i11::g.43161::m.43161  0.523839
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733            0.469819
# TRINITY_DN1272_c0_g1::TRINITY_DN1272_c0_g1_i7::g.6644::m.6644       0.168933
# stat
# <numeric>
#   TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747              5.72724
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.744::m.744              5.46313
# TRINITY_DN1275_c0_g1::TRINITY_DN1275_c0_g1_i9::g.6875::m.6875       -4.79519
# TRINITY_DN17372_c0_g1::TRINITY_DN17372_c0_g1_i11::g.43161::m.43161   4.77881
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733             4.69065
# TRINITY_DN1272_c0_g1::TRINITY_DN1272_c0_g1_i7::g.6644::m.6644       -4.63343
# pvalue
# <numeric>
#   TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747            1.02079e-08
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.744::m.744            4.67812e-08
# TRINITY_DN1275_c0_g1::TRINITY_DN1275_c0_g1_i9::g.6875::m.6875      1.62521e-06
# TRINITY_DN17372_c0_g1::TRINITY_DN17372_c0_g1_i11::g.43161::m.43161 1.76335e-06
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733           2.72339e-06
# TRINITY_DN1272_c0_g1::TRINITY_DN1272_c0_g1_i7::g.6644::m.6644      3.59661e-06
# padj
# <numeric>
#   TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747            0.000173545
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.744::m.744            0.000397663
# TRINITY_DN1275_c0_g1::TRINITY_DN1275_c0_g1_i9::g.6875::m.6875      0.007494691
# TRINITY_DN17372_c0_g1::TRINITY_DN17372_c0_g1_i11::g.43161::m.43161 0.007494691
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733           0.009260071
# TRINITY_DN1272_c0_g1::TRINITY_DN1272_c0_g1_i7::g.6644::m.6644      0.010191005

# overall, there are not many consistent genes between generations. 

head(resAM_OW)  

# summary statistics for ambient vs ocean warming and acidification

summary(resAM_OW)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ plotting individual genes


### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)

d <-plotCounts(dds, gene="TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733", intgroup = (c("treatment","generation")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7)

p

# the above graph is showing how warming is driving the difference in expression! you can see the higher count numbers for the OW and OWA treatment groups in comparison to AM and OA. Note that this graph is looking at the expression of one specific gene of choice (as indicated in line 350)


####################################################### looking between generations

#################### MODEL NUMBER 2 - subset to focus on effect of treatment for each generation

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)

dim(dds)
# [1] 67916    38

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_sub <- subset(dds, select = generation == 'F0')
dim(dds_sub)

# Perform DESeq2 analysis on the subset
dds_sub <- DESeq(dds_sub)
# [1] 20598 [2] 12

resultsNames(dds_sub)

res_F0_OWvAM <- results(dds_sub, name="treatment_OW_vs_AM", alpha=0.05)

res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

summary(res_F0_OWvAM)


# out of 20598 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3006, 15% #up-regulated
# LFC < 0 (down)     : 1835, 8.9%
# outliers [1]       : 16, 0.078%
# low counts [2]     : 0, 0%
# (mean count < 26)


### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds_sub, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment","generation")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

# this is showing differential expression for generation F0 ^ 

# sideways volcano plot (MA plot)
plotMA(res_F0_OWvAM, ylim=c(-4,4))


######################################## HEATMAP 

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

# By environment
vsd <- vst(dds_sub, blind=FALSE)

topgenes <- head(rownames(res_F0_OWvAM),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_sub)[,c("generation","treatment")])
pheatmap(mat, annotation_col=df)
pheatmap(mat, annotation_col=df, cluster_cols = F)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VENN DIAGRAM~~~~~~~~~~~~~

#################################################################

#### PLOT OVERLAPPING DEGS IN VENN EULER DIAGRAM

#################################################################

# For OW vs AM
res_F0_OWvAM <- results(dds_sub, name="treatment_OW_vs_AM", alpha=0.05)
res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

summary(res_F0_OWvAM)
res_F0_OWvAM <- res_F0_OWvAM[!is.na(res_F0_OWvAM$padj),]
degs_F0_OWvAM <- row.names(res_F0_OWvAM[res_F0_OWvAM$padj < 0.05,])

# For OA vs AM
res_F0_OAvAM <- results(dds_sub, name="treatment_OA_vs_AM", alpha=0.05)
res_F0_OAvAM <- res_F0_OAvAM[order(res_F0_OAvAM$padj),]
head(res_F0_OAvAM)

summary(res_F0_OAvAM)
res_F0_OAvAM <- res_F0_OAvAM[!is.na(res_F0_OAvAM$padj),]
degs_F0_OAvAM <- row.names(res_F0_OAvAM[res_F0_OAvAM$padj < 0.05,])

# For OWA vs AM
res_F0_OWAvAM <- results(dds_sub, name="treatment_OWA_vs_AM", alpha=0.05)
res_F0_OWAvAM <- res_F0_OWAvAM[order(res_F0_OWAvAM$padj),]
head(res_F0_OWAvAM)

summary(res_F0_OWAvAM)
res_F0_OWAvAM <- res_F0_OWAvAM[!is.na(res_F0_OWAvAM$padj),]
degs_F0_OWAvAM <- row.names(res_F0_OWAvAM[res_F0_OWAvAM$padj < 0.05,])

library(eulerr)

# Total
length(degs_F0_OAvAM)  # 520
length(degs_F0_OWvAM)  # 4841 
length(degs_F0_OWAvAM)  # 3742

# Intersections
length(intersect(degs_F0_OAvAM,degs_F0_OWvAM))  # 387
length(intersect(degs_F0_OAvAM,degs_F0_OWAvAM))  # 340
length(intersect(degs_F0_OWAvAM,degs_F0_OWvAM))  # 2585

intWA <- intersect(degs_F0_OAvAM,degs_F0_OWvAM)
length(intersect(degs_F0_OWAvAM,intWA)) # 308

# Number unique

520-387-340+308 # 101 OA
4841-387-2585+308 # 2177 OW 
3742-340-2585+308 # 1125 OWA

387-308 # 79 OA & OW
340-308 # 32 OA & OWA
2585-308 # 2277 OWA & OW


# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("OA" = 101, "OW" = 2177, "OWA" = 1125, "OA&OW" = 79, "OA&OWA" = 32, "OW&OWA" = 2277, "OA&OW&OWA" = 308))


plot(fit1,  lty = 1:3, quantities = TRUE)
# lty changes the lines

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))


#cross check
2177+2277+308+79 # 4841, total OW
1125+2277+308+32 # 3742, total OWA
101+32+79+308    # 520, total OA






