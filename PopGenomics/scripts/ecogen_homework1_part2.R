# ecogen homework 1: part 2: GEA Analysis

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("")

list.files()

### read in selection statistics (these are chi^2 distributed)

s<-npyLoad("allRS_poly_hw1.selection.npy")

as.data.frame(s)

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,2))
names(pval) = c("p_PC1","p_PC2")

## read positions
p <- read.table("allRS_poly_hw1_mafs.sites",sep="\t",header=T, stringsAsFactors=T)
dim(p)

p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)

# get all the outliers with p-values below some cutoff
cutoff=1e-3   

outliers_PC2 <- p_filtered[which(pval$p_PC2<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff?
dim(outliers_PC2)[1]


# write them out to a file
write.table(outliers_PC2,
            "allRS_poly_outliers_PC2.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F,
            )




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

COV <- as.matrix(read.table("allRS_poly_hw1.cov"))

PCA <- eigen(COV)

as.data.frame(PCA$vectors)

data=as.data.frame(PCA$vectors)
data=data[,c(2:3)] # the second number here is the number of PC axes you want to keep

write.table(data,
            "allRS_poly_genPC2_3.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F,
            )
#eol="/n"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# You made need to use "install.packages" if you don't have some of the below libraries already

library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)

setwd("")

bio <- getData("worldclim",var="bio",res=12)

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)

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