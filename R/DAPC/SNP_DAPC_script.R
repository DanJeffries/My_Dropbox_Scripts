#### DAPC analysis ####

## 6676 SNPs, pure crucian. Analysis in '/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/populations/crucian_only_M8_N8/'

library(adegenet)

## CWD

setwd('/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_V3/populations/crucian_only_M8_N8/DAPC')


body <- read.genepop('./batch_1.gen') ## read in file

## Cluster the data with PCA
grp <- find.clusters (body, max.n.clust = 50)

## Chose number of PCs to retain: 100
## Chose 4 clusters - BIC scores looked interesting. Definitely support for 4 cluster, but also perhaps 18ish clusters
## which is about the number of populations.

table.value(table(pop(body), grp$grp), col.lab=paste("inf", 1:8), row.lab=paste("Pop", body$pop.names)) ## take a look at the cluster membership table

## pretty conclusive - PRO and V belong to their own, separate clusters, and the rest fall into 2 clusters.

## Do quick test DAPC and perform spline analysis to find optimal number of PCs to retain.
dapc.test <- dapc(body, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(dapc.test)
## 6 is the optimal number suggested

##So now do the real thing, chosing 6 PCs
dapc.real <- dapc(body, grp$grp)
## The eigenvalues for the linear discriminants point to using just one LD. So try this.

## Take a look at the results

scatter(dapc.real, scree.da=T, posi.da="topleft", bg = "gray80", pch= 20, solid=.6, cex=3)
## the plot is interesting, although doesn't tell me much about the which colour is for which group. no idea how to get this


########### second DAPC with 2 LDs chosen ##############

dapc2LDs <- dapc(body, grp$grp)

scatter(dapc2LDs, scree.da=T, posi.da="bottomright", bg = "gray80", pch= 20, solid=.6, cex=3)

### I think this is a much better representation - on the single LD plot, group 1 looks too close to groups 2 and 3.
## Groups are 1 = V (River Danube), 2 = UK,DEN,S.SWEDEN, 3 = N.SWEDEN, FINLAND, POLAND, 4 = PRO (River Don)



###### DAPC with N. Europe subset ###########

## Because of the divergence between N. Europe and PRO and V populations, I need to re-run with these taken out. 
## However, taking these out is not as easy as it sounds - I think I need to re-run populations with no V or PRO.. then redo


body$ind.names
N.EU <- read.genepop('../Pure_cruc_N_EU_r_0.8_P15/batch_1.gen')



grp <- find.clusters (N.EU, max.n.clust = 50)
## Chose 100
## BIC scores suggest 10 clusters

table.value(table(pop(N.EU), grp$grp), col.lab=paste("inf", 1:10), row.lab=paste("Pop", N.EU$pop.names)) ## take a look at the cluster membership table



dapc.test <- dapc(N.EU, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(dapc.test)
## choose 6 PCs according to a-score


dapcNEU <- dapc(N.EU, grp$grp)
## 6 PCs, 3 LDs
scatter(dapcNEU, scree.da=T, posi.da="bottomright", bg = "gray80", pch= 20, solid=.6, cex=3)
