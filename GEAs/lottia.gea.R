####################### Lottia GEAs #######################

#Genotype-environment association tests using WGS data from Lottia gigantea
#
# inspired by this tutorial- https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html
#
# Adapted and/or written by ES Nielsen, UC Davis, USA
#
# Code is provided as is, without support 

###########################################################

LIB <- c("rgbif", "ggplot2", "gridExtra", "knitr", "raster", 
         "ade4", "rworldmap", "maptools", "rasterVis", "rgdal","rgeos", "sdmpredictors")
#for(i in LIB) { install.packages(i, repos="http://ftp.sun.ac.za/ftp/pub/mirrors/cran.za.r/") ; library(i, character.only=T) }
for(i in LIB) { library(i, character.only=T) }

library(usdm)
library(psych)    
library(vegan)
library(adegenet)
library(dplyr)
library(fmsb)
library(gsl)
library(raster)
library(sdmpredictors)
library(adegenet)
library(qvalue)
library(tidyverse)

######################## DOWNLOAD PREDICTOR VARIABLES

#download landscape features
#bio1= Annual Mean Temperature
#bio5= Max Temperature of Warmest Month
#bio6= Min Temperature of Coldest Month
#bio7= Air Temperature Range
#bio12 = Annual Mean Precipitation
a.contp <- load_layers( layercodes = c("WC_bio5", "WC_bio6", "WC_bio1", "WC_bio7", "WC_bio12") , equalarea=FALSE, rasterstack=TRUE)

#download marine features
o.contp <- load_layers( layercodes = c("BO22_temprange_ss", "BO22_tempmean_ss", "BO22_salinitymean_ss", "BO21_salinityrange_ss", "BO22_ppmean_ss", "BO22_dissoxmean_ss", "BO22_ph", "BO22_cloudmean", "BO22_chlomean_ss", "BO22_curvelmean_ss", "BO22_damean") , equalarea=FALSE, rasterstack=TRUE)

#resample because air and ocean vars have different resolution/extent
a<-extent(-180, 180, -90, 90) #input layer extent
o<-extent(-180, 180, -90, 90)
extent_list<-list(a, o)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)
reso<-res(a.contp) #choose layer you want to keep resolution
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=a.contp@crs) #choose layer crs you want to keep
a.c.r.2 <-resample(a.contp, s2, method="ngb") #resample by s2
o.c.r.2 <-resample(o.contp, s2, method="ngb") #resample by s2
contemp.r.2=stack(a.c.r.2, o.c.r.2) #stack resampled layers

#Crop layers to West Coast extent
WC.ext <- extent(-130, -107, 20, 50)
env.curr.c <- crop(contemp.r.2, WC.ext)

# Load sample site lat longs
lat_longs <- read_excel("lat.longs.xlsx", col_types = c("text", "numeric", "numeric"))

# Extract env values for our sample locations:

envdata <- raster::extract(x=env.curr.c, y=cbind(lat_longs$x, lat_longs$y))
data<- cbind2 (x=lat_longs, y=envdata)

#Extract when there's NAs
env.curr.c <- stack(env.curr.c)
sample_raster_NA <- function(r, xy){
  apply(X = xy, MARGIN = 1, 
        FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
  
}
test<-lapply(env.curr.c@layers, function(a_layer) sample_raster_NA(a_layer, cbind(lat_longs$x, lat_longs$y)))
envdat.whole <-as.data.frame(do.call(cbind, test))
data<- cbind2 (x=lat_longs, y=envdat.whole)

pairs.panels(envdata, scale=T)
env.cors.all <- cor(envdat.whole, method = c("spearman"))
env.cors.all <- as.data.frame(env.cors.all)
write.table(env.cors.all, "spearman.corrs.all.env.txt")


vif(env.curr.c)

VARSEL <- c("BO22_damean", "BO22_tempmean_ss", "BO22_ph")
contemp.mar.3vars <- stack(subset(env.curr.c, VARSEL))
vif(contemp.mar.3vars)

Variables      VIF
1      BO22_damean 1.513805
2 BO22_tempmean_ss 2.241261
3          BO22_ph 2.456132

test<-lapply(contemp.mar.3vars@layers, function(a_layer) sample_raster_NA(a_layer, cbind(lat_longs$x, lat_longs$y)))
envdat.3 <-as.data.frame(do.call(cbind, test))
colnames(envdat.3) <- c("BO22_damean","BO22_tempmean_ss","BO22_ph")

writeRaster(contemp.mar.3vars, "GEA.env.3vars.tif", format="GTiff", overwrite=TRUE)
write.table(envdat.3, "GEA.env.per.site.txt")

######################## RUN GEA WITH LFMM

str_name<-'GEA.env.3vars.tif' 
envdat<-imported_raster=raster(str_name)

# convert vcf to lfmm
setwd("~/Desktop/PD_stuffies/GEAs/LFMM")

library(lfmm)
library(LEA)

### Get gen data
gen.input <- read.lfmm("by_pop_0.05_pctind0.5.lfmm")

#check number of samples and snps with dimensions
dim(gen.input)

#check for missing genotypes
sum(is.na(gen.input))

### Get env data
lfmm.env <- read.table("GEA.env.per.site.txt", header=TRUE)

str(lfmm.env) 

env.1 <-(lfmm.env[,1]) #damean

#run lfmm
lottia.lfmm <- lfmm_ridge(gen.input, lfmm.env, K=3) ## change K as you see fit

#calculate test statistics for the predictors
lottia.pv <- lfmm_test(Y=gen.imp, X=lfmm.env, lfmm=lottia.lfmm, calibrate="gif")

# this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values
names(lottia.pv)

# Let’s look at the genomic inflation factor (GIF)
# An appropriately calibrated set of tests will have a GIF of around 1
lottia.pv$gif

#NOTE: Changing the value of K influences the GIF, so additional tests using the “best” value of K +/- 1 may be needed in some cases. See Francois et al. (2016) for more details.

# plot unadjusted & adjusted p-values
# you want to see a relatively flat histogram (most loci not under selection) 
# with a peak near zero, indicating candidate adaptive markers

hist(lottia.pv$pvalue[,1], main="Unadjusted p-values")        
hist(lottia.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

# Let's change the GIF and readjust the p-values:
zscore <- lottia.pv$score[,1]   # zscores for first predictor, we only have one in our case...
(gif <- lottia.pv$gif[1])       ## d.fault GIF for the first predictor

new.gif1 <- 0.9  #choose gif as you see fit

# Manual adjustment of the p-values:
adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)

#plot manually adj p-values
hist(adj.pv1, main="REadjusted p-values (GIF=1.0)") #this gives us a better peak at 0,
                                                    #meaning not too conservative

#lottia.pv.adj <- cbind2 (x=lottia.pv, y=adj.pv1)

#convert the adjusted p-values to q-values

#with lfmm chosen GIF
lottia.qv <- qvalue(lottia.pv$calibrated.pvalue)$qvalues

#with our manually adjusted GIF
lottia.qv <- qvalue(adj.pv1)$qvalues

lottia.qv.df <-as.data.frame(lottia.qv)

length(which(lottia.qv < 0.1)) ## how many SNPs have an FDR < 10%?

(lottia.FDR.1 <- colnames(gen.imp)[which(lottia.qv < 0.1)]) ## identify which SNPs these are

# re-do with adjsted p-vals for each env variable
# (basically re-run p-val script with df=env var number)

#run lfmm older version
sites_environ_matrix<-as.matrix(env.input[,-c(1:5)])
write.env(sites_environ_matrix, "sites_environ_matrix.env")

test_lfmm <- lfmm("pruned.lfmm", "sites_environ_matrix.env", K=2, rep=1, project = "new")



######################## RUN GEA WITH RDA
library(vegan)

setwd("~/Desktop/PD_stuffies/GEAs/RDAs")

lottia_geno <- read.table("ALL_2.rda.geno.txt", header=TRUE)
dim(lottia_geno)
#[1] 703925     19

lottia_freq<-t(lottia_geno)

lottia_env <- read.table("~/Desktop/PD_stuffies/GEAs/RDAs/NEW/env.xy.NEW.txt", header=TRUE)
lottia_env <- read.table("GEA.env.per.site.txt", header=TRUE)
lottia_env <- env.xy.NEW[,2:4]

#Run the RDA
lottia.rda <- rda(lottia_freq ~ ., data=lottia_env, scale=T)

#Look at R2, or variance explained. i.e. a adjusted R^2 of 0.0002 
# means that our constrained ordination explains about 0.02% of the variation.
RsquareAdj(lottia.rda)

#$r.squared
$r.squared
[1] 0.2281082

$adj.r.squared
[1] 0.07372982

#on adjusted R2: Our constrained ordination explains about 7% of the variation; this low explanatory power is not surprising given that we expect that most of the SNPs in our dataset will not show a relationship with the environmental predictors (e.g., most SNPs will be neutral)

R.sum <- summary(lottia.rda)
R.sum$cont   # Prints the "Importance of components" table
R.sum$cont$importance[2, "RDA1"]
#  0.119951
R.sum$cont$importance[2, "RDA2"]
#  0.05654108

# observe the eigenvalues for the constrained axes
# reflect the variance explained by each canonical axis

summary(eigenvals(lottia.rda, model = "constrained"))
screeplot(lottia.rda)

# test for sign. environ. variables
anova(lottia.rda, by="terms", permu=200)
anova(lottia.rda, by="axis", perm.max=500) # test axes for significance
anova(lottia.rda, perm=999) #test for overal significance

#check collinearity 
vif.cca(lottia.rda)


# do forward selection
fwd.sel <- ordiR2step(rda(lottia_freq ~ 1, data = lottia_env, scale=T), # lower model limit (simple!)
                      scope = formula(lottia.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) # change to TRUE to see the selection process!

fwd.sel$call

#plot RDA 
#SNPs are red at the center of the plot
#population are black 
#vectors are env predictors
pdf("lottia.rda.NEW.12.pdf")
plot(lottia.rda, scaling=3)
text(lottia.rda, display="sites", col=1, scaling=3)
dev.off()

pdf("lottia.rda.NEW.13.pdf")
plot(lottia.rda, choices=c(1,3), scaling=3)#axis 1 and 3
text(lottia.rda, display="sites", col=1, scaling=3)
dev.off()

### plot RDA color by population
library(vegan)
library(cowplot)

scrs <- scores(lottia.rda, display = c("species", "sites", "bp"), choices = 1:4, scaling = 2)
scrs

snp_centroids <- data.frame(scrs$species)
snp_centroids

snp_centroids$PC_names <- rownames(snp_centroids) 

aem_continuous_arrows <- data.frame(scrs$biplot)
aem_continuous_arrows

rownames(aem_continuous_arrows) <- c("DA","SST","PH")

aem_continuous_arrows$aem_number <- rownames(aem_continuous_arrows)
aem_continuous_arrows

baseplot<-plot(lottia.rda, scaling = 2)

mult <- attributes(baseplot$biplot)$arrow.mul

sitenames <- sites$Abb

### identify candidates

#get first 3 loadings (stored as species in the RDA object)
load.rda <- summary(lottia.rda)$species[,1:3]

#look at dist (should be close to normally dist) - those on tails likely candidates
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

cand1 <- outliers(load.rda[,1], 3) #3 for 3 std dev 
cand2 <- outliers(load.rda[,2], 3) 
cand3 <- outliers(load.rda[,3], 3) 

ncand<-length(cand1) + length (cand2) + length (cand3)
ncand
[1] 1232

#lottia.rda.cand <- c(names(cand1), names(cand2), names(cand3)) ## j.st the names of the candidates

#Next, we'll organize our results by making one data frame with the axis, SNP name, loading, & correlation with each predictor:
cand1.df <- cbind.data.frame(rep(1, times = length(cand1)), names(cand1), unname(cand1)); colnames(cand1.df) <- c("axis", "snp", "loading")
cand2.df <- cbind.data.frame(rep(2, times = length(cand2)), names(cand2), unname(cand2)); colnames(cand2.df) <- c("axis", "snp", "loading")
cand3.df <- cbind.data.frame(rep(3, times = length(cand3)), names(cand3), unname(cand3)); colnames(cand3.df) <- c("axis", "snp", "loading")

cand <- rbind(cand1.df, cand2.df, cand3.df)
cand$snp <- as.character(cand$snp)

cand.mat <- matrix(nrow=(ncand), ncol=3)  # ncol = number of predictors
colnames(cand.mat) <- c("DA", "SST", "pH")
                        
#Let's add in the correlations of each candidate SNP with the  environmental predictors:
n<-dim(lottia_env)[2]
foo <- matrix(nrow=ncand, ncol=3)  # ncol= number of env predictors
colnames(foo) <- colnames(lottia_env)

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- lottia_freq[,nam]
  foo[i,] <- apply(lottia_env,2,function(x) cor(x,snp.gen))
}

#table of candidates snp with loading on each axis and correlation with env predictors
cand <- cbind.data.frame(cand,foo)  
head(cand)

#investigate the candidates
#are there snp associated with several axis? if yes remove them
n_dupli<-length(cand$snp[duplicated(cand$snp)])
n_dupli
if (n_dupli>=1){cand<-cand[!duplicated(cand$snp)]}

#Next, we'll see which of the predictors each candidate SNP is most strongly correlated with:
n<-dim(cand)[2]
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,n+1] <- names(which.max(abs(bar[3:n]))) # gives the variable
  cand[i,n+2] <- max(abs(bar[3:n]))              # gives the correlation
}

colnames(cand)[n+1] <- "predictor"
colnames(cand)[n+2] <- "correlation"
head(cand)
table(cand$predictor) 
write.table(cand, "lottia.rda.NEW.outs.txt", quote=FALSE, sep=" ", row.names=FALSE)

#Let's look at RDA plots again, but this time focus in on the SNPs in the ordination space. 
#We'll color code the SNPs based on the predictor variable that they are most strongly correlated with. 

sel <- cand$snp
env <- cand$predictor
variables<-levels(as.factor(cand$predictor))
for (i in 1: length(variables))
{env[env==variables[i]] <- i} #associates a color with each env predictor

# color by predictor:
col.pred <- rownames(lottia.rda$CCA$v) # pull the SNP names and put them in a vector of colour

for (i in 1:length(sel)) {           
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("LOTGIsca",col.pred)] <- '#f1eef6' # non-candidate SNPs make them transparent
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c(1: length(variables))

# axes 1 & 2
pdf("lottia.rda.NEW.outs.1_2.pdf")
plot(lottia.rda, type="n", scaling=1, xlim=c(-1,1), ylim=c(-1,1))
points(lottia.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=1)
points(lottia.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=1)
text(lottia.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomleft", legend=variables, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
dev.off()


# axes 1 & 3
pdf("lottia.rda.NEW.outs.1_3.pdf")
plot(lottia.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(lottia.rda, display="species", pch=21, cex=1.5, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(lottia.rda, display="species", pch=21, cex=1.5, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(lottia.rda, scaling=3, display="sites", col="black", cex=1.5)
text(lottia.rda, scaling=3, display="bp", col="#0868ac", cex=1.5)
legend("bottomleft", legend=c("DA","SST","pH"), bty="n", col="gray32", pch=21, cex=1.5, pt.bg=bg)
dev.off()


#check environmental predictors are most strongly correlated with the first three RDA axes
intersetcor(lottia.rda)[,1:3]

################ Running GEO RDA 
lottia_xy <- env.xy.NEW[,5:6]

#Run the RDA
lottia.xy.rda <- rda(lottia_freq ~ ., data=lottia_xy, scale=T)

#Look at R2, or variance explained. i.e. a adjusted R^2 of 0.0002 
# means that our constrained ordination explains about 0.02% of the variation.
RsquareAdj(lottia.xy)

#$r.squared
$r.squared

$adj.r.squared



############### Running pcRDAs (to account for pop structure)
gen.pca <- rda(lottia_geno, scale=T)
screeplot(gen.pca, main = "Screeplot: Eigenvalues of Genetic Variables")

loadings <- scores (gen.pca, display = 'species', scaling = 0)
#copied the first PC scores and added them as column to env var txt input

lottia_env_pc <- read.table("env.vars.pc.txt", header=TRUE)
pairs.panels(lottia_env_pc, scale=T) #shows that PC is not correlated with any of our env vars, so we can use in RDA

lottia.rda.pc <- rda(lottia_freq ~ BO22_damean + BO22_tempmean_ss + BO22_ph + Condition(PC1), data=lottia_env_pc, scale=T)

RsquareAdj(lottia.rda.pc)
$r.squared
#[1] 0.1774947

$adj.r.squared
#[1] 0.02381132

anova(lottia.CA.rda, perm=999)
anova(lottia.rda, by="terms", permu=999)


############### Running dbRDAs (to account for IBD)
library(adespatial)
library("ejanalysis/proxistat")

lottia_env_xy <- read.table("env.vars.coords.txt", header=TRUE)

Coor=lottia_env_xy[,4:5]

plot(Coor, asp=1)

## Create dbMEM
nsites.per.group <- c(19)
dbmem <- create.dbMEM.model(coord=Coor, nsites=nsites.per.group)
summary(dbmem)

#Combine dbMEM and environmental variable as an unique dataset named Y.
Y=cbind(dbmem, lottia_env_xy[,1:3])

#look for correlations between MEMs & env vars
cor(Y) 
pairs.panels(Y, scale=T) # dbmem1 is correlated with pH and SST

#remove MEM1 as it is correlated with other vars
db.pred <-Y[,2:6]

#Run the RDA
lottia.rda.db <- rda(lottia_freq ~ ., data=db.pred, scale=T)

RsquareAdj(lottia.rda.db)
#$r.squared
#[1] 0.320163

#$adj.r.squared
#[1] 0.05868726


######### RDA with both PC to account for pop struc & MEMs for IBD
lottia_env_pc <- read.table("env.vars.pc.txt", header=TRUE)
lottia_pred_all <-cbind(lottia_env_pc, db.pred)
#selecting vars to use (not keeping mems 1&2 due to collinearity)
rda.pred.vif <-lottia_pred_all %>% select("BO22_damean", "BO22_tempmean_ss", "BO22_ph", "PC1", "dbMEM.3")

pairs.panels(rda_pred_all, scale=T) #no colinearity issues 

lottia.rda.all <- rda(lottia_freq ~  dbMEM.3 + BO22_damean + BO22_tempmean_ss + BO22_ph + Condition(PC1), data=rda.pred.vif, scale=T)

RsquareAdj(lottia.rda.all)
$r.squared
[1] 0.2256661

$adj.r.squared
[1] 0.01896889

vif.cca(lottia.rda.all)
PC1          dbMEM.3      BO22_damean BO22_tempmean_ss 
2.215047         1.634503         2.010476         3.042247 
BO22_ph 
2.622553

###################### Run RDA with just CA sites & protection level

#read new env data
env.xy.prot.CA <- read.delim("~/Desktop/PD_stuffies/GEAs/RDAs/env.xy.prot.CA.txt")
CA.env<-env.xy.prot.CA[c(8,2:4)]

#create new gen data
dim(lottia_geno)
#[1] 703925     19

CA.gen<-lottia_geno[,1:15]
dim(CA.gen)

CA.freq<-t(CA.gen)

#Run the RDA
lottia.CA.rda <- rda(CA.freq ~ ., data=CA.env, scale=T)

#Look at R2, or variance explained. i.e. a adjusted R^2 of 0.0002 
# means that our constrained ordination explains about 0.02% of the variation.
RsquareAdj(lottia.CA.rda)

#$r.squared
0.314373

#$adj.r.squared
0.04012217

#plot RDA 
#SNPs are red at the center of the plot
#population are black 
#vectors are env predictors
pdf("CA.rda.12.pdf")
plot(lottia.CA.rda, scaling=3)
text(lottia.CA.rda, display="sites", col=1, scaling=3)
dev.off()

pdf("CA.rda.13.pdf")
plot(lottia.CA.rda, choices=c(1,3), scaling=3)#axis 1 and 3
text(lottia.CA.rda, display="sites", col=1, scaling=3)
dev.off()

######################## COMPARE LFMM & RDA

intersect(lottia.FDR.1, lottia.rda.cand) ## f.und by both LFMM and RDA

setdiff(lottia.FDR.1, lottia.rda.cand)   # unique to LFMM


# to get snp list with sequential numbering
setwd("~/Desktop/PD_stuffies/GEAs/baypass")
snp.list <- read.table("snp.list.txt", header = F)

snp.numbs <- snp.list %>% mutate(observation = 1:n()) %>% select(observation, everything())
names(snp.numbs)[1] <- "number"

#read baypass fst out list
bpass.cov.list <- read.table("bpass.cov.out.list.txt", header = T)
names(bpass.cov.list)[1] <- "number"
names(bpass.cov.list)[2] <- "snp"

#read baypass env out lists
bpass.list <- read.table("bpass.NEW.env.out.list.txt", header = T)
names(bpass.list)[2] <- "number"

#read lfmm out lists
setwd("~/Desktop/PD_stuffies/GEAs/LFMM")
lfmm.list <- read.table("lfmm.uniq.NEW.outs.txt", header = F)
names(lfmm.list)[1] <- "number"

#change column names
lfmm.snp.numbs <- merge(snp.numbs, lfmm.list, by  = "number") 
bpass.snp.numbs <- merge(snp.numbs, bpass.list, by  = "number") 

#read rda out lists
setwd("~/Desktop/PD_stuffies/GEAs/RDAs")
rda.list <- read.table("lottia.rda.NEW.outs.txt", header = T)

#change column names
names(bpass.snp.numbs)[2] <- "snp"
names(lfmm.snp.numbs)[2] <- "snp"

#get overlapping snps for GEAs
lfmm.bpass.overlap <- merge(lfmm.snp.numbs, bpass.snp.numbs, by  = "number") 
lfmm.rda.overlap <- merge(lfmm.snp.numbs, rda.list, by  = "snp") 
bpass.rda.overlap <- merge(bpass.snp.numbs, rda.list, by  = "snp")

#get overlapping SNPs with fst version of Baypass
lfmm.bp.overlap <- merge(lfmm.snp.numbs, bpass.cov.list, by  = "number") 
rda.bp.overlap <- merge(bpass.cov.list, rda.list, by  = "snp") 
bpass.bp.overlap <- merge(bpass.snp.numbs, bpass.cov.list, by  = "snp")

setwd("~/Desktop/PD_stuffies/GEAs")
write.table(lfmm.bpass.overlap, "NEW.lfmm.bpass.overlap.txt", quote=FALSE, sep=" ", row.names=FALSE)

write.table(lfmm.rda.overlap, "NEW.lfmm.rda.overlap.txt", quote=FALSE, sep=" ", row.names=FALSE)

write.table(bpass.rda.overlap, "NEW.bpass.rda.overlap.txt", quote=FALSE, sep=" ", row.names=FALSE)

write.table(lfmm.bp.overlap, "NEW.lfmm.bp_cov.overlap.txt", quote=FALSE, sep=" ", row.names=FALSE)

write.table(rda.bp.overlap, "NEW.rda.bp_cov.overlap.txt", quote=FALSE, sep=" ", row.names=FALSE)

write.table(bpass.bp.overlap, "NEW.bpass.bp_cov.overlap.txt", quote=FALSE, sep=" ", row.names=FALSE)


#get total snp list to remove from neutral dataset
bpass<-bpass.snp.numbs[2]
lfmm<-lfmm.snp.numbs[2]
rda<-rda.list[2]
bp<-bpass.cov.list[2]

full.out.list<-rbind(bpass, bp, lfmm, rda)

write.table(full.out.list, "NEW.full.out.list.txt", quote=FALSE, sep=" ", row.names=FALSE)

#import sites list (from angsd, which was edited to be one column with chrom_pos)
sites <- read.table("sites_ALL_2_maf0.05_pctind0.5_maxdepth15_R2.pruned", header = F)

#remove outs from sites list
sites.neut <- sites[!(sites$V1 %in% full.out.list$snp),]

write.table(sites.neut, "sites_ALL_2_R2.pruned.neut.txt", quote=FALSE, sep=" ", row.names=FALSE)

#creat new regions file for neutral only
#import regions list
regs <- read.table("regions_ALL_2_maf0.05_pctind0.5_maxdepth15", header = F)

#import nuet sites file
sites <- read.table("sites_ALL_2_R2.pruned.neut.NEW.txt", header = T)

#subset
reg.neut <- subset(regs, V1 %in% sites$V1)

#write
write.table(reg.neut, "regions_ALL_2_R2.pruned.neut.NEW.txt", quote=FALSE, sep=" ", row.names=FALSE)


#get outlier snp list (take overlap snps; chosen by > 2 methods)
outs1<-as.data.frame(lfmm.bpass.overlap[,2])
outs2<-as.data.frame(lfmm.rda.overlap[,1])
outs3<-as.data.frame(bpass.rda.overlap[,1])

outs4<-as.data.frame(lfmm.bp.overlap[,2])
outs5<-as.data.frame(rda.bp.overlap[,1])
outs6<-as.data.frame(bpass.bp.overlap[,1])

names(outs1)[1] <- "snp"
names(outs2)[1] <- "snp"
names(outs3)[1] <- "snp"
names(outs4)[1] <- "snp"
names(outs5)[1] <- "snp"
names(outs6)[1] <- "snp"

#merge dataframes of overlapping snps
library(tidyverse)

#put all data frames into list
df_list <- list(outs1, outs2, outs3, outs4, outs5, outs6)

#merge all data frames in list
outs_df <- df_list %>% reduce(full_join, by='snp')


#remove duplicate snps
outs.final<-outs_df[!duplicated(outs_df$snp), ]
outs.final<-as.data.frame(outs.final)

write.table(outs.final, "NEW.doub.out.list.txt", quote=FALSE, sep=" ", row.names=FALSE)

## R script to get overlapping SNPs btwn GEA.all GEA.CA.n.c

setwd("~/Desktop/PD_stuffies/GEAs")
all.list <- read.table("NEW.doub.out.list.txt", header = T)

setwd("~/Desktop/PD_stuffies/R_work/CA.n.c.outs")
CA.nc.list <- read.table("CA.n.c.360.list.txt", header = T)

setwd("~/Desktop/PD_stuffies/R_work/harv.outs")
harv.list <- read.table("BP_C2.prot.harv.out.list.txt", header = T)


ALL.CA.nc.list <- intersect(all.list, CA.nc.list) ## f.und by both LFMM and RDA

library(dplyr)
ALL.CA.nc.list <-inner_join(all.list, CA.nc.list)