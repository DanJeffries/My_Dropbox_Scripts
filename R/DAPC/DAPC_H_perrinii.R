library(adegenet)
library(reshape2)
library(ggplot2)
library(scatterplot3d)


## RUNNING THE DAPC ANALYSIS ==================================================================

setwd('/home/djeffrie/Data/RADseq/H_perrinii/Populations_kept_largepops_p12_r0.8_m6/')

body <- read.genepop('./batch_1.gen') ## read in genepop SNP file (outputted by Stacks populations)


### I ran steps 1 - 4 below twice, once for K = 2 and once for k = 3

# 1. Cluster the data with PCA =======================================================================================
grp <- find.clusters (body, max.n.clust = 50)
#(retained all 120 plotted PCs in this step)

table.value(table(pop(body), grp$grp), col.lab=paste("inf", 1:11), row.lab=paste("Pop", popNames(body))) ## take a look at the cluster membership table

# 2. Do quick test DAPC and perform spline analysis to find optimal number of PCs to retain. =========================
dapc.test <- dapc(body, n.da=100, n.pca=20, grp$grp)
temp <- optim.a.score(dapc.test)

## 6 PCs recommended

# 3. Now do the real DAPC, chosing 6 PCs =============================================================================
dapc.real <- dapc(body, grp$grp)
## The eigenvalues for the linear discriminants point to using just one LD. So try this.

# 4. PLOTTING ========================================================================================================

cols = c("blue", "green") #, "green")

### Scatter / Density  

scatter(dapc.real, col = cols, scree.da=T, posi.da="topleft", bg = FALSE, pch= 20, solid=.6, cex=3)

### Plotting structure like plot 

compoplot(dapc.real, col = cols)


### Write the posterior probabilities of assignment to each K for individuals

write.table(dapc.real$posterior, file = "./K2_assignment_posteriors.txt",
    sep = "\t", quote = F)

dev.off()


### Plot an additional normal PCA

body_noNA <- tab(body, NA.method="mean")

pca1 <- dudi.pca(body_noNA,scannf=FALSE,scale=FALSE)

# get pop identities

body_pop <- as.character(body$pop)

pops = vapply(strsplit(body_pop,"_"), `[`, 1, FUN.VALUE=character(1))

pca_info = data.frame(sample_ID = indNames(body), pop = pops,  PC1 = pca1$li$Axis1, PC2 = pca1$li$Axis2, K = dapc.real$assign)


palette(rainbow(3))

plot(pca_info$PC1, pca_info$PC2, col = pca_info$K, pch = as.character(pca_info$pop))





