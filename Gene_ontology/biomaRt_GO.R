### R script to do GO

setwd("~/Desktop/PD_stuffies/R_work/outs_GO")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("biomaRt")

library("biomaRt")
library(topGO)

listMarts()

# look at top 10 databases      
head(biomaRt::listMarts(host = "https://www.ensembl.org"), 10)  

#collect gene names from biomart

mart <- biomaRt::useMart(biomart = "metazoa_mart",
                         dataset = "lgigantea_eg_gene",
                         host = 'metazoa.ensembl.org')


# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id",
                                        "go_id"), mart = mart)

#examine result
head (GTOGO)

#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]

# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))

#examine result
head (geneID2GO)

#select genes of interest from the gene list
library(dplyr)

#read LDannot gene id list (I created this by edited the LD-annot output in excel and just keeping column of gene IDs)
LDannot_gene_ids <- read_excel("LDannot.gene.ids.xlsx")
int.genes.df <-as.data.frame(LDannot_gene_ids)

#keep only unique gene IDs from biomaRt
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))

#filter all gene list to only include those from LD-annot
#cands.df <-subset(all.gene.df, all.genes %in% LDannot_gene_ids$gene_id)

cands.chr <-as.character(int.genes.df[[1]])

# compare the chr lists of cands & all genes and turn into factor as to whether they overlap or not
int.genes <- factor(as.integer(all.genes %in% cands.chr))
names(int.genes) = all.genes

# #create GO object (for running in topGO!)
library(topGO)

# GO for Biological Process
go.obj.bp <- new("topGOdata", ontology='BP'
              , allGenes = int.genes # this is factor containing interesting genes (gene IDs for outliers)
              , annot = annFUN.gene2GO # we use this function as we supply our own annotations 
              , gene2GO = geneID2GO # this is our object we just created with biomaRt
)
resultFisher <- getSigGroups(go.obj.bp, test.stat)
resultWeight <- getSigGroups(go.obj.bp, test.stat2)
allRes <- GenTable(go.obj.bp, classic = resultFisher, 
                   weight = resultWeight, orderBy = "weight", 
                   ranksOf = "classic",topNodes=length(resultFisher@score),numChar=100)
filtRes.BP <- allRes[allRes$classic<0.05 & allRes$Significant>2,]
filtRes.BP$Ontology <- "BP"

# GO for Molecular Function
go.obj.mf <- new("topGOdata", ontology='MF'
              , allGenes = int.genes # this is factor containing interesting genes (gene IDs for outliers)
              , annot = annFUN.gene2GO # we use this function as we supply our own annotations 
              , gene2GO = geneID2GO # this is our object we just created with biomaRt
)
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(go.obj.mf, test.stat)
test.stat2 <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(go.obj.mf, test.stat2)
allRes <- GenTable(go.obj.mf, classic = resultFisher, 
                   weight = resultWeight, orderBy = "weight", 
                   ranksOf = "classic",topNodes=length(resultFisher@score),numChar=100)
filtRes.MF <- allRes[allRes$classic<0.05 & allRes$Significant>2,]
filtRes.MF$Ontology <- "MF"

# GO for Cellular Component
go.obj.cc <- new("topGOdata", ontology='CC'
                 , allGenes = int.genes # this is factor containing interesting genes (gene IDs for outliers)
                 , annot = annFUN.gene2GO # we use this function as we supply our own annotations 
                 , gene2GO = geneID2GO # this is our object we just created with biomaRt
)
resultFisher <- getSigGroups(go.obj.cc, test.stat)
resultWeight <- getSigGroups(go.obj.cc, test.stat2)
allRes <- GenTable(go.obj.cc, classic = resultFisher, 
                   weight = resultWeight, orderBy = "weight", 
                   ranksOf = "classic", topNodes=length(resultFisher@score))
filtRes.CC <- allRes[allRes$classic<0.05 & allRes$Significant>2,]
filtRes.CC$Ontology <- "CC"

##Output all three
filt.all <- rbind(filtRes.BP,filtRes.MF,filtRes.CC)
write.csv(filt.all,"GO.new.GEA.outs.csv")



### OLD- not used!! 
#run topGO to get molecular function of outlier genes (using Fisher test)
results <- runTest(go.obj, algorithm = "elim", statistic = "fisher")


#create a table with gene ontology
results.tab <- GenTable(object = go.obj, elimFisher = results)

#write as CSV
write.csv(results.tab,"GO.output.csv")

