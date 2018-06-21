library(adegenet)
library(reshape2)
library(ggplot2)
library(scatterplot3d)


## RUNNING THE DAPC ANALYSIS ==================================================================

setwd('/home/djeffrie/Data/RADseq/Hori/Stacks_new/Populations_all_kept_final/')

body <- read.genepop('./batch_1.gen') ## read in genepop SNP file (outputted by Stacks populations)

# 1. Cluster the data with PCA =======================================================================================
grp <- find.clusters (body, max.n.clust = 50)

#(retained all PCs in this step)

table.value(table(pop(body), grp$grp), col.lab=paste("inf", 1:11), row.lab=paste("Pop", popNames(body))) ## take a look at the cluster membership table

# 2. Do quick test DAPC and perform spline analysis to find optimal number of PCs to retain. =========================
dapc.test <- dapc(body, n.da=100, n.pca=50, grp$grp)
temp <- optim.a.score(dapc.test)

## 6 PCs recommended

# 3. Now do the real DAPC, chosing 6 PCs =============================================================================
dapc.real <- dapc(body, grp$grp)
## The eigenvalues for the linear discriminants point to using just one LD. So try this.

# 4. Plotting ========================================================================================================

## Take a look at the results in a scatter
scatter(dapc.real, col = rainbow(10), scree.da=T, posi.da="topleft", bg = "gray90", pch= 20, solid=.6, cex=3)

## Trying 3D plot

coordinates <- as.data.frame(dapc.real$ind.coord)

mycols = rainbow(12)

scat3d <- scatterplot3d(coordinates$LD1, coordinates$LD3, coordinates$LD2, 
                        color = "black", 
                        bg = mycols[dapc.real$assign], 
                        pch = 21, 
                        cex.symbols = 2, 
                        xlab = "DA1", 
                        ylab = "DA3",
                        zlab = "DA2",
                        angle = 70, 
                        type = "h",
                        scale.y = 1, 
                        lwd= 0.5)

## Not that useful. 


### Plotting structure like plot ================================================================

## Adegenet has the "compoplot" function which will do this automatically. 

svg(filename = "../../Population_structure/Structure_like_plot_v2.svg",
    width = 15,
    height = 5)
compoplot(dapc.real, col = rainbow(10))
dev.off()


