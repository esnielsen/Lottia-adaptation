## Script to run PCA on pcangsd output

library(dplyr)
library(ade4)
library(ggplot2)
library(viridis)

library(tidyverse)
library(dplyr)
library(cluster)
library(factoextra)

 setwd("~/Desktop/PD_stuffies/R_work/PCAs")

# Pop colors
 
CA.c.n.cols=c("#436EEE", "#6959CD","#009ACD","#B9D3EE","#40E0D0","#458B00","#66CD00", "#00FF00","#8B7500","#FFFF00", "#FFD700","#EE7600", "#CD3700")
 
CA.cols=c("#436EEE", "#6959CD","#009ACD","#B9D3EE","#40E0D0","#458B00","#66CD00", "#00FF00","#8B7500","#FFFF00", "#FFD700","#EE7600", "#CD3700", "#FF1493","#FF83FA")

#ALL.cols=c("#68228B", "#9F79EE", "#CD6090", "#8B1C62", "#FF83FA", "#FF1493", "#CD3700", "#FFA500", "#FFD700", "#FFFF00", "#EEEE00", "#00FF00", "#66CD00", "#458B00", "#40E0D0","#B9D3EE", "#009ACD", "#6959CD", "#436EEE")

ALL.cols.NEW=c("#436EEE", "#6959CD","#009ACD","#B9D3EE","#40E0D0","#458B00","#66CD00", "#00FF00","#8B7500","#FFFF00", "#FFD700","#EE7600", "#CD3700", "#FF1493","#FF83FA", "#8B1C62", "#CD6090", "#9F79EE", "#68228B")

#harvesting factor
siteharv <- c("N", "Y", "N", "Y", "N", "Y", "Y", "N", "N", "Y", "N", "N", "Y")
  
## Download pca file and info file from cluster
pca.mat<-read.table(paste0("ALL2_OUTS_CAnc_maf0.05_pctind0.5_maxdepth15_pruned.cov.pca"), header=T)

meta_info<-read.table("all.filt.info.txt", header=T) # !!! here make sure that the first column is not read as row names
meta_info<-read.table("CA.nc.info.txt", header=T) 
lat<-meta_info[,3]
#harv<-meta_info[,8]
pop<-meta_info[,2]

pop=factor(pop, levels=c("KR", "FR", "BB", "DB", "HP", "CR", "SB", "VN", "GP", "CA", "WA", "SH", "SC", "CB", "PB", "SR", "BT", "BA"), ordered = TRUE)


#### Plot and color by group
# NB: you can also change the group to column 2 above and it plots by individual

#pca of all sites
pca<-ggplot(pca.mat,aes(x=PC1,y=PC2)) + geom_point(aes(color=lat)) + scale_colour_manual(values = rev(ALL.cols.NEW))+theme_minimal()

#pca of CA sites
pca<-ggplot(pca.mat,aes(x=PC1,y=PC2)) + geom_point(aes(color=lat)) + scale_colour_manual(values = rev(CA.cols))+theme_minimal()

#pca of CA sites w/protection shape
pca<-ggplot(pca.mat,aes(x=PC1,y=PC2)) + geom_point(aes(color=lat, shape=harv), size=2.2) + scale_colour_manual(values = rev(CA.cols))+theme_minimal()


pdf("ALL_2.pca.pruned.pdf")
pca
dev.off()

#### Plot with clusters
pca_transform = as.data.frame(pca.mat[,1:2])

fviz_nbclust(pca_transform, kmeans, method = 'wss')
k = 3
kmeans_pca = kmeans(pca_transform, centers = k, nstart = 50)

fviz_cluster(kmeans_pca, data = pca.mat,
             palette = ALL.cols.NEW, 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

### plot with color by latitutde

# Get latitude
lat <- (meta_info[,3])

ggplot(pca.mat,aes(x=PC1,y=PC2)) + geom_point(aes(color=lat)) + scale_color_viridis(discrete = TRUE, option='plasma')

### plot with plotly
plotly::plot_ly(pca.mat, type = "scatter", mode = "markers", x=pca.mat$PC1, y=pca.mat$PC2, text=bam_x_info$id)

#Get variance explained
vars_transformed <- apply(pca.mat, 2, var)
vars_transformed/sum(vars_transformed)


#### Plot with UMAP
library(umap)
library(M3C)

cov<-as.matrix(read.table("ALL_2_maf0.05_pctind0.5_maxdepth15_pruned.cov"))
umap(cov,labels=as.factor(x_info$pop),controlscale=TRUE,scale=3, colvec = ALL.cols)


#### plot as in Marine Genomics
cov<-as.matrix(read.table("ALL_2_maf0.05_pctind0.5_maxdepth15.cov"))
#e<-eigen(cov)
#plot(e$vectors[,1:2])

# this gives % variance explained by each PC
#e$values/sum(e$values)

#### plot the .covMAt and .ibsMat from angsd pca script #
#MDS
name <- "ALL_angsd.pca_maf0.05_pctind0.5_maxdepth15.ibsMat"
m <- as.matrix(read.table(name))
mds <- cmdscale(as.dist(m))
plot(mds,lwd=2,ylab="Dist",xlab="Dist",main="multidimensional scaling")

#PCA
name <- "ALL.2_angsd.pca_maf0.05_pctind0.5_maxdepth15.covMat"
m <- as.matrix(read.table(name))
e <- eigen(m)
plot(e$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 2",main="Principal components",col=rep(1:19,each=30),pch=16)


PC<-as.data.frame(e$vectors[,1:2])
info <- cbind(PC, x_info)
ggplot(info,aes(x=V1,y=V2)) + geom_point(aes(color=lat)) + scale_color_viridis(discrete = TRUE, option='plasma')

#### neighbor joining tree (doesn't really work)
 heatmap / clustering / trees
name <- "angsdput.ibsMat" # or covMat
m <- as.matrix(read.table(name))


plot(ape::nj(m))
plot(hclust(dist(m), "ave"))

#just in case order is different, use rownames of the pca (coming from bam.list) to order the info matrix
# bam_name<-as.data.frame(rownames(pca.mat))
# colnames(bam_name)<-colnames(all.filt.info)[1]
# bam_x_info<-left_join(bam_name,all.filt.info)

#Get group names 
#NB: make sure the pca.mat rows are in same order as column 1 of bam info!!!
# group<-as.factor(bam_x_info[,3])


#p<-ggplot(pca.mat,aes(x=PC1,y=PC2)) + geom_point(aes(color=lat)) + scale_color_viridis(discrete = TRUE, option='plasma') +theme_minimal()

#pp<-p+geom_point()

#pp


