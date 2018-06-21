## From adagenet tutorial ##

install.packages("adegenet")

library(adegenet)
library(ape)

adegenetTutorial("genomics"):

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_all/populations_r07_p_30_m8/")

## made the PLINK input files using vcftools from the batch_1.vcf file outputted by stacks

alldata <- read.PLINK('batch_1_plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)


##check data

pop(alldata) ## contains individual names instead - change this:

## note - this is different from alldata$pop


### Adding my own population names ####

pops <- read.delim('pops.txt', header = F)
samples <- read.delim("sample_names.txt", header = F)

## made this file using cut from the pop_codes.txt file used in "populations". 

## check everything
pop(alldata) <- pops$V1
indNames(alldata) <- samples$V1
pop(alldata) ## looks good
nInd(alldata)
indNames(alldata)
other(alldata) ## nothing in here
ploidy(alldata)
head(alleles(alldata)) ## nothing in here

## all looks good!

head(as.matrix(alldata)) ## can look at the data as a matrix, and it can be subsetted the same way using [indx.row, idx.col] where rows are individuals and columns are loci



plot(alldata, posi = 'topleft', cex - 0.5) ## can't get pop names on here yet - the yaxis is annoyingly unsuppressable. 
axis(2, at = seq(1,170,3), labels = rev(pop(alldata)[seq(1,length(pop(alldata)), 3)]), las = 2, cex = 0.01) ## complicated axis plotting - need to reverse the populations, because the function plots down to up. And also need to only plot every 3rd pop name because the axis gets too cramped
?plot


### look at the number of alleles per locus

myfreq <- glMean(alldata) ## looking at one allele only here

hist(myfreq, proba = T, col="gold", xlab = 'Allele Freqs', main = " Distribution of second allele frequencies")

## to look at alleles frequencies:

myfreq <- c(myfreq, 1-myfreq)
hist(myfreq, proba = T, col="gold", xlab = 'Allele Freqs', main = " Distribution of allele frequencies", nclass = 20)
temp <-density(myfreq, bw = .05)
lines(temp$x, temp$y*2, lwd = 3)


## Subset data

pop(alldata) <- pops$V1

x <- seploc(alldata, n.block = 10, parallel=F)

y <- seppop(alldata)

crudata <- repool(y$BF, y$BOR, y$CAKE, y$CALK, y$COP, y$MOAT, y$OBY, y$OU, y$PED, y$POLEN, y$PRO, y$RM, y$SD, y$SK, y$STYV, y$STEC, y$TROM, y$TU, y$V, y$WEN)
names(x)
names(y)

## Make tree based on genetic distance

DISTS <- lapply(x, function(e) dist(as.matrix(e)))
names(DISTS)
D <- Reduce("+", DISTS)

par(mar = c(1,1,1,1))
plot(nj(D), cex =0.5)

unrootedtree <- nj(D) ## unroot

rooted <- (root(unrootedtree, interactive = T)) ## or root

plot(rooted, cex = 0.5)

pca1 <- glPca(alldata)

scatter(pca1, posi = "bottomleft")

### PCA ########

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0,col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2,posi='bottomright', inset=.05,ratio=.3)

##### Neighbour Joining tree #####

tree <- nj(D)

## plot the tree
plot(tree, cex = 0.5)

## plot tree with the PCA results on it.
plot(tree, typ='fan', cex = 0.2)
tiplabels(pch=20,col=myCol, cex=2)



