###################################
#  Selection scans for red spruce #
###################################

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("~/Documents/Github/Ecological_Genomics/Fall_2023/pcangsd/")

list.files()

### read in selection statistics (these are chi^2 distributed)

s<-npyLoad("allRS_poly.selection.npy")

head(s)

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,1))
head(pval)
names(pval) = c("p_PC1", "p_PC2")
head(pval)

## read positions
p <- read.table("allRS_poly_mafs.sites",sep="\t",header=T, stringsAsFactors=T)

dim(p)
head(p)

# filter data frame 

p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)

# How many sites got filtered out when testing for selection? Why?

## make manhattan plot
plot(-log10(pval$p_PC1),
     col=p_filtered$chromo,
     xlab="Position",
     ylab="-log10(p-value)",
     main="Selection outliers: pcANGSD e=1 (K2)")

# We can zoom in if there's something interesting near a position...

plot(-log10(pval$p_PC1[2e05:2.01e05]),
     col=p_filtered$chromo, 
     xlab="Position", 
     ylab="-log10(p-value)", 
     main="Selection outliers: pcANGSD e=1 (K2)")

# get the contig with the lowest p-value for selection
sel_contig <- p_filtered[which(pval==min(pval$p_PC1)),c("chromo","position")]
sel_contig

# get all the outliers with p-values below some cutoff
cutoff=1e-4   # equals a 1 in 5,000 probability
outlier_contigs <- p_filtered[which(pval<cutoff),c("chromo","position")]
outlier_contigs

outlier_contig <- outlier_contigs[which(outlier_contigs$position>0),]

# how many outlier loci < the cutoff?
dim(outlier_contig)[1]

# how many unique contigs harbor outlier loci?
length(unique(outlier_contig$chromo))

write.table(unique(outlier_contig$chromo),
            "allRS_poly_PC1_outlier_contigs.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ october 4

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("")

list.files()

### read in selection statistics (these are chi^2 distributed)

s<-npyLoad("allRS_poly.selection.npy")

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,1))
names(pval) = c("p_PC1","p_PC2")

## read positions
p <- read.table("allRS_poly_mafs.sites",sep="\t",header=T, stringsAsFactors=T)
dim(p)

p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)

# get all the outliers with p-values below some cutoff
cutoff=1e-3   

outliers_PC1 <- p_filtered[which(pval$p_PC1<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff?
dim(outliers_PC1)[1]


# write them out to a file
write.table(outliers_PC1,
            "allRS_poly_outliers_PC1.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F)

# read in: 

COV <- as.matrix(read.table("allRS_poly.cov"))

PCA <- eigen(COV)

data=as.data.frame(PCA$vectors)
head(data)

data=data[,c(1:2)] # the second number here is the number of PC axes you want to keep

# write out scores using write table function:
write.table(data,
            "allRS_poly_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

### Getting Climate data: BIOCLIM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)


bio <- getData("worldclim",var="bio",res=10)

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)

head(coords)

#The chunk below refers to your bamlist file that you transferred during last week's PCA/admixture analysis.  It should be the same one you want to use here -- if your sample list for analysis changes in the future, you'll need a different bamlist!
  
names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

angsd_coords <- merge(pops, coords, by.x="Ind", by.y="Tree")

points <- SpatialPoints(angsd_coords[c("Longitude","Latitude")])

clim <- extract(bio,points)

angsd_coords_clim <- cbind.data.frame(angsd_coords,clim)
str(angsd_coords_clim)

#angsd_coords_clim has the bioclim variables we are interested in as well as 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA Analysis with climate data 
# Make the climate PCA:

clim_PCA = PCA(angsd_coords_clim[,15:33], graph=T) #column 15 is where bioclim 1 begins 

# if a population is located near the error for bio1, for example, it comes from a region of high seasonal variation

# Get a screeplot of climate PCA eigenvalues

fviz_eig(clim_PCA)

# What is the climate PCA space our red spruce pops occupy?

fviz_pca_biplot(clim_PCA, 
                geom.ind="point",
                col.ind = angsd_coords_clim$Latitude, 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                title="Climate PCA (Bioclim)",
                legend.title="Latitude")

# Which variables show the strongest correlation on the first 2 climate PC axes?

dimdesc(clim_PCA)[1:2]

# ^ bio12 has the highest PC correlation, at 0.937


# __________________________________________________________________ exploring most significant bioclim variables

# Replace "XX" with your bio variable most significant on climate PC1:

write.table(scale(angsd_coords_clim["bio12"]),
            "allRS_bio12.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)


# Replace "YY" with your bio variable most significant on climate PC2:  

write.table(scale(angsd_coords_clim["bio10"]),
            "allRS_bio10.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)

# at this point, we export these tables to the server to run our GEA analysis. 










