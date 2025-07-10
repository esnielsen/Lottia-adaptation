# Lottia-adaptation

Pipeline used to analyze Lottia gigantea whole genome sequencing data, for the publication titled "The effects of latitudinal gradients, climatic anomalies, and size-selective harvesting on the adaptive potential of an intertidal gastropod" by Nielsen, Erica; Walkes, Samuel; Sones, Jacqueline; Fenberg, Phillip; Paz Garcia, David; Grosberg, Rick; Sanford, Eric; Bay, Rachael.

Corresponding miro board of full pipeline can be found here: https://miro.com/welcomeonboard/THE0eXpPUXMxcVFnNHd6OENtYXdMZk1NRTd6SEFzbHdFZTl2a25OSHh6TmQ1OFF5NFcyQU9lUkpRck96YlZObXwzNDU4NzY0NTI0MDc1NzU3Mjc0fDI=?share_link_id=714947409411

The following is a combination of bash, R, and python code written/adapted by Erica Nielsen with Rachael Bay at UC Davis. We ran the following on the Pittsburgh Supercomputing Center's Bridges2 servers.



## Part 1 - Download data, call SNPs, & generate population diversity metrics

As this work is an extension of the publication titled "Pushed waves, trailing edges, and extreme events: Eco-evolutionary dynamics of a geographic range shift in the owl limpet, Lottia gigantea" (which looks at neutral genomic variation in the species), the beginning steps of the pipeline are outlined in this repository: https://github.com/esnielsen/Lottia_range_expansion. Follow all steps in Part 1 and Part 2 of the repository to get SNPs and genetic diversity indices used within this repository. 



## Part 2 - Genotype-environment association analyses (GEAs)
 
Within the GEAs folder of this repo, you will find an R code called *lottia.gea.R* which is used to do the following:
- downlod environmental predictor variables from Bio-Oracle and WorldClim
- Run GEA with LFMM
- Run GEA with RDAs
- Get subset of outlier SNPs which overlap between at least 2 outlier tests

There is also a script called *bpass.sim.pod.R* which is used with the following bash scripts:
-*baypass1*
-*baypass2*

We also ran IBE tests, which was done with the *ibd.ide.R* R code. 

To generate PAI and SGV and map both of these metrics, use the following code: *plot.svg.pai.R*


## Part 3 - Leading-edge outliers

Within the this folder of the repo, you will find the following scripts:
-*baypass C2* = uses the BayPass C2 model to identify outliers (this same script is also used to identify outliers by comparing harvesting groups)
-*baypass C2* = uses Population Branch statistic to identify outliers


## Part 4 - Harvesting genomic differentiation

In the harvesting folder you will find the scripts:
-*pop.het.R* = generate expected heterozygosity values
-*cmh.prot.tests.R* = run CMH test to identify outliers


 ## Part 5 - Gene ontology

 In the GOterms folder you will find scripts that identify the functional roles of outliers detected above
 *LD annot scripts*
-*biomaRt.GO.R* = gets GO terms from LDannot output


  ## Part 6 - HGAMs

This folder contains data and scripts to run the HGAMs. 

 - *LD annot scripts* = raw data from MARINe
 - *LD annot scripts* = wrangled data as input for GAMs
 - *LD annot scripts* = R code to run GAMs
