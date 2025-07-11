#### Script to run cmh test for allele diffs between protected/non-protected areas

setwd("~/Desktop/PD_stuffies/R_work/protect.outs")
library(poolSeq)
library(qvalue)

#read baypass input file
bpass.in<- read.table("~/Desktop/PD_stuffies/GEAs/baypass/by_pop_0.05_pctind0.5_maxdepth15.mafs.pruned.baypass.txt", header=F)

#transform (rows=pops, cols=inds)
cmh.in <- t(bpass.in)
dim(cmh.in)

#get minor allele counts for protected area sites
p0<- cmh.in[c(1,5,9,17,23,29),]

#get major allele counts for protected area sites
P0<- cmh.in[c(2,6,10,18,24,30),]

#get minor allele counts for non-protected area sites
np0<- cmh.in[c(3,7,11,19,25,27),]

#get major allele counts for non-protected area sites
NP0<- cmh.in[c(4,8,12,20,26,28),]

# perform CMH-test for all empirical loci
p.values <- cmh.test(A0=P0, a0=p0, At=NP0, at=np0, min.cov=1, max.cov=0.99, min.cnt=1, log=TRUE)


# plot -log10 transformed p-values
plot(p.values, main="CMH Manhattan plot", xlab="SNP", ylab="-log10(p)", pch=".")

all_qvals <- qvalue(p.values)
FDR <- 0.1
which(all_qvals$qvalues < 0.1)

[1] 342056= LOTGIsca_21_3312170

pval_cutoff <-max(all_qvals$pvalues[all_qvals$qvalues <= FDR])
1.050089e-07 

