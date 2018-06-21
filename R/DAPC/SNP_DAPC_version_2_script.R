#### DAPC analysis ####

######## WHOLE EU DATASET ##############

## "Dataset contains only first SNP of each locus. Whole dataset contained 26142 SNPs
##- 13344 SNPs produced using the --write-single-SNP flag
##- Mean number of SNPs per tag = 1.96.", pure crucian. 
## Analysis in '/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_analysis_V2'

library(adegenet)


## All outputs saved here: '/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_analysis_V2/populations_r_0.7_P_16_m_8_single_SNP/DAPC'

setwd('/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_analysis_V2/populations_all_pops_blacklisted/')


body <- read.genepop('./batch_1.gen') ## read in genepop SNP file (outputted by Stacks populations)

## Cluster the data with PCA
grp <- find.clusters (body, max.n.clust = 50)
## Chose number of PCs to retain: 150
## Chose 10 clusters based on BIC scores


table.value(table(pop(body), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("Pop", body$pop.names)) ## take a look at the cluster membership table

## Do quick test DAPC and perform spline analysis to find optimal number of PCs to retain.
dapc.test <- dapc(body, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(dapc.test)
## 7 is the optimal number suggested by a-score

##So now do the real thing, chosing 7 PCs
dapc.real <- dapc(body, grp$grp)
## The eigenvalues for the linear discriminants point to using just one LD. So try this.

## Take a look at the results

scatter(dapc.real, scree.da=T, posi.da="topleft", bg = "gray80", pch= 20, solid=.6, cex=3)
## the plot is interesting, although doesn't tell me much about the which colour is for which group. no idea how to get this


########## 5 cluster coercion ###############


## Cluster the data with PCA
grp <- find.clusters (body, max.n.clust = 50)

## Chose number of PCs to retain: 200
## Chose 5 clusters based on BIC scores


table.value(table(pop(body), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("Pop", body$pop.names)) ## take a look at the cluster membership table

## Do quick test DAPC and perform spline analysis to find optimal number of PCs to retain.
dapc.test <- dapc(body, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(dapc.test)
## 7 is the optimal number suggested by a-score

##So now do the real thing, chosing 7 PCs
dapc.real <- dapc(body, grp$grp)
pops <- body$pop.names
par(mar = c(5.1,4.1,5,1.1))
complot<- compoplot(dapc.real, lab ="", las = 3, cex.lab = 0.5, posi = list(x=0, y=1.2), cleg =0.7)
axis(1, at = complot, labels = body$ind.names, cex = 0.5)
?compoplot

str(body$tab)

## Take a look at the results

scatter(dapc.real, scree.da=T, posi.da="topleft", bg = "gray80", pch= 20, solid=.6, cex=3)
## the plot is interesting, although doesn't tell me much about the which colour is for which group. no idea how to get this



########### second DAPC with 2 LDs chosen ##############

dapc2LDs <- dapc(body, grp$grp)

scatter(dapc2LDs, scree.da=T, posi.da="bottomright", bg = "gray80", pch= 20, solid=.6, cex=3)

### I think this is a much better representation - on the single LD plot, group 1 looks too close to groups 2 and 3.
## Groups are 1 = V (River Danube), 2 = UK,DEN,S.SWEDEN, 3 = N.SWEDEN, FINLAND, POLAND, 4 = PRO (River Don)



############ THE ABOVE IS NOT REPRODUCIBLE - EACH TIME I RUN IT THROUGH I GET A DIFFERENT
############ NUMBER OF CLUSTERS SUGGESTED BY THE BIC SCORES. HOWEVER, WHEN I CONSTRAIN TO 5 CLUSTERS
############ THE RESULTS ARE A BIT MORE STABLE BUT THE CLUSTER SHARING IS STILL VERY PLASTIC BETWEEN RUNS!



##### Northern Europe Only ###########

### This is a separate Stacks analysis - with a new catalog, sstacks and population run performed. See working dir below path for outputs
## SNPs were conditional on being present in 70% or over in at least 16 out of the 18 populations used 


setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_NEU_only/NEU_r_0.7_p_16_m8_single_SNP/")


## Dataset contains only first SNP of each locus. Whole dataset contained 28988 loci - only 7614 of these conained 1 or more SNPs
## So 7614 SNPs produced using the --write-single-SNP flag, out of 14043 SNPs across all loci.
## Note, Mean number of SNPs per tag = 1.845.

### all outputs in ('/media//dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_NEU_only/NEU_r_0.7_p_16_m8_single_SNP/DAPC_outputs/')


N.EU <- read.genepop('./batch_1.gen')

N.EU$pop.names

grp <- find.clusters (N.EU, max.n.clust = 50)
## Chose 100
## Chose 7 clusters
## BIC scores support several numbers of clusters

table.value(table(pop(N.EU), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("Pop", N.EU$pop.names)) ## take a look at the cluster membership table



dapc.test <- dapc(N.EU, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(dapc.test)



dapcNEU <- dapc(N.EU, grp$grp)
## choose 7 PCs according to a-score
## Chose 2 LDs
scatter(dapcNEU, scree.da=T, posi.da="bottomright", bg = "gray80", pch= 20, solid=.6, cex=3)


########## ALTHOUGH THE NUMBER OF CLUSTERS IS MORE STABLE HERE, THE RESULTS ARE STILL VARIABLE 
########## ACROSS RUNS.... WHY IS THIS? 


#-----------------------------------------------------------------------------------------------------------
####### Checking DIYABC 1000 SNP dataset ###############



## All outputs saved here: '/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/Pure_crucian_analysis_V2/populations_r_0.7_P_16_m_8_single_SNP/DAPC'

setwd('~/Dropbox/PhD/Dans_PhD_Shared/Papers//Phylogeography paper//Final_data_files/')


body <- read.fstat('Pure_cru_1000_rndm_FSTAT.dat') ## read in genepop SNP file (outputted by Stacks populations)

pops <- read.delim("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Data/RAD/DIYABC/pops.txt", header = F)
pops
body$pop.names
pop(body) <- pops$V1

## Cluster the data with PCA
grp <- find.clusters (body, max.n.clust = 50)
## Chose number of PCs to retain: 150
## Chose 10 clusters based on BIC scores


table.value(table(pop(body), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("Pop", body$pop.names)) ## take a look at the cluster membership table

## Do quick test DAPC and perform spline analysis to find optimal number of PCs to retain.
dapc.test <- dapc(body, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(dapc.test)
## 7 is the optimal number suggested by a-score

##So now do the real thing, chosing 7 PCs
dapc.real <- dapc(body, grp$grp)
## The eigenvalues for the linear discriminants point to using just one LD. So try this.

## Take a look at the results

scatter(dapc.real, scree.da=T, posi.da="topleft", bg = "gray80", pch= 20, solid=.6, cex=3)
## the plot is interesting, although doesn't tell me much about the which colour is for which group. no idea how to get this



