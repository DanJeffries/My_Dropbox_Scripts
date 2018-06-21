library(adegenet)
library(calibrate)
library(ggplot2)

#-----------------------------------------------------------------------------------------------------
#######################################################################################################
###################################### All_RAD_data ##################################################
######################################################################################################
#-----------------------------------------------------------------------------------------------------

setwd("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter 4 Hybridisation and introgression/data/RAD")

#### Doing PCA in adegenet #####

RADdata <- read.PLINK('plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T) # read in data

## Add population info
RAD_pops <- read.delim("PCS_pops.txt", header = F)
pop(RADdata) <- RAD_pops$V1

## Add sample info
RAD_samples <- read.delim("sample_names.txt", header = F)
indNames(RADdata) <- RAD_samples$V1

pop(RADdata)

pca1 <- glPca(RADdata, nf = 10) ## do PCA with max 10 retained PCs
barplot(pca1$eig[1:20],main="PCA eigenvalues", col=heat.colors(20)) ## PCA eigens

## Get PCA plotting data for ggplot2
RAD_PCA_coords<- pca1$scores
RAD_PCA_dat <- as.data.frame(RAD_PCA_coords) ## make it a dataframe
RAD_PCA_dat$pops <- pop(RADdata) ## add pop info

RAD_SNP_loadings_mat <- pca1$loadings ## add loc names
row.names(RAD_SNP_loadings_mat) <- RADdata$loc.names ## am I sure these are in the right order??
RAD_SNP_loadings <- as.data.frame(RAD_SNP_loadings_mat)

## Get PC percentage variances
perc_var_PC1 = (pca1$eig[1]/sum(pca1$eig))*100
perc_var_PC2 = (pca1$eig[2]/sum(pca1$eig))*100
perc_var_PC3 = (pca1$eig[3]/sum(pca1$eig))*100
perc_var_PC4 = (pca1$eig[4]/sum(pca1$eig))*100
perc_var_PC5 = (pca1$eig[5]/sum(pca1$eig))*100

paste("PC1 (",round(perc_var_PC1, 2),"%)")


#------------------------------------------------------------------------------------------------
                                 ###### PLOTTING THE PCA #########


### The below PCAs are very nice - and seem to pull out the axes of variation well ##

pca <- qplot(PC1,PC2, data = RAD_PCA_dat, color = pops, cex = 2, xlab = paste("PC1 (",round(perc_var_PC1, 2),"%)"), ylab = paste("PC2 (",round(perc_var_PC2, 2),"%)"))
pca + annotate("text", x = RAD_PCA_dat$PC1, y = RAD_PCA_dat$PC2, label = row.names(RAD_PCA_dat), cex = 4)
## PC1 - Carassius - Cyrpinus
## PC2 - C.carassius -> C.auratus -> C. gibelio


#### Getting allele loadings PC1####
orderedPC1<- RAD_SNP_loadings[order((RAD_SNP_loadings$Axis1)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot(((orderedPC1$Axis1)^2)[c(1:500)], names.arg = row.names(orderedPC1[c(1:500),]), cex.names = 0.3, las = 3)



#### Getting allele loadings PC2 ####
orderedPC2<- RAD_SNP_loadings[order((RAD_SNP_loadings$Axis2)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot(((orderedPC2$Axis2)^2)[c(1:100)], names.arg = row.names(orderedPC2[c(1:100),]), cex.names = 0.5, las = 2)


## PC2 vs PC3

pca <- qplot(PC2,PC3, data = RAD_PCA_dat, color = pops, cex = 2, xlab = paste("PC2 (",round(perc_var_PC2, 2),"%)"), ylab = paste("PC3 (",round(perc_var_PC3, 2),"%)"))
pca + annotate("text", x = RAD_PCA_dat$PC2, y = RAD_PCA_dat$PC3, label = row.names(RAD_PCA_dat), cex = 4)

#### Getting allele loadings PC3 ####
orderedPC3<- RAD_SNP_loadings[order((RAD_SNP_loadings$Axis3)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot(((orderedPC3$Axis3)^2)[c(1:100)], names.arg = row.names(orderedPC3[c(1:100),]), cex.names = 0.5, las = 2)

## PC3 - C.carassius Hap1 - C.carassius Hap2


## PC3 Vs PC4 

pca <- qplot(PC3,PC4, data = RAD_PCA_dat, color = pops, cex = 2, xlab = paste("PC3 (",round(perc_var_PC3, 2),"%)"), ylab = paste("PC4 (",round(perc_var_PC4, 2),"%)"))
pca + annotate("text", x = RAD_PCA_dat$PC3, y = RAD_PCA_dat$PC4, label = row.names(RAD_PCA_dat), cex = 4)

#### Getting allele loadings PC4 ####
orderedPC4 <- RAD_SNP_loadings[order((RAD_SNP_loadings$Axis4)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot(((orderedPC4$Axis4)^2)[c(1:100)], names.arg = row.names(orderedPC4[c(1:100),]), cex.names = 0.5, las = 2)

## PC4 - C.auratus - C.gibelio


## PC4 Vs PC5 
pca <- qplot(PC4,PC5, data = RAD_PCA_dat, color = pops, cex = 2, xlab = paste("PC4 (",round(perc_var_PC4, 2),"%)"), ylab = paste("PC5 (",round(perc_var_PC5, 2),"%)"))
pca + annotate("text", x = RAD_PCA_dat$PC4, y = RAD_PCA_dat$PC5, label = row.names(RAD_PCA_dat), cex = 4)

#### Getting allele loadings PC3####
orderedPC5<- RAD_SNP_loadings[order((RAD_SNP_loadings$Axis5)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot(((orderedPC5$Axis5)^2)[c(1:100)], names.arg = row.names(orderedPC5[c(1:100),]), cex.names = 0.5, las = 2)

## PC5 - Looks like its within C. gibelio and is less informative

#### Finally - output the snp loadings to a file for manuscript ####

top100PC1 <- row.names(orderedPC1[c(1:100),])
top100PC2 <- row.names(orderedPC2[c(1:100),])
top100PC3 <- row.names(orderedPC3[c(1:100),])
top100PC4 <- row.names(orderedPC4[c(1:100),])
top100PC5 <- row.names(orderedPC5[c(1:100),])

top100s <- cbind(top100PC1,top100PC2,top100PC3,top100PC4,top100PC5) ## combine these into a matrix to be written

write.matrix(top100s, "top100_snp_loadings.csv", sep = '\t') ## write to file



#-----------------------------------------------------------------------------------------------------
#######################################################################################################
################################## All_microsat_data ##################################################
######################################################################################################
#-----------------------------------------------------------------------------------------------------

setwd("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter 4 Hybridisation and introgression/data/Microsats/")

### Doing PCA in Adegenet, following tutorial on page 51 ###

Microsats_all <- read.fstat("Car_refined_final.dat")  ## get data and make genind obj
Micro_names <- read.delim("NAMES_final.txt", header = F) ## get names 
indNames(Microsats_all) <- Micro_names$V1 ## add names to genind

X <- scaleGen(Microsats_all, missing="mean") ## Scale gen replaces na's with best guess data - required for pca analysis 

pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=6) ## Doing PCA
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50)) ## Plotting the PC eigenvalues

Micro_PCA<- pca1$li ## These are the cordinates for each individual on each PC
allele_loadings <- pca1$c1 ## These are the allele loadings for every allele on every PC

#-------------------------------------------------------------------------------------------------------
                         ########## Plotting the PCA ##############

## PC1 Vs PC2


pca <- qplot(Axis1,Axis2, data = Micro_PCA, color = pop(Microsats_all), cex = 2)
pca + annotate("text", x = Micro_PCA$Axis1, y = Micro_PCA$Axis2, label = row.names(Micro_PCA), cex = 3)

## PC1 (Axis1) Pulls out Carassius carassius - carassius goldfish variation
## PC2 Pulls out Carassius - cyprinus variation


#### Getting allele loadings PC1####
orderedPC1<- allele_loadings[order((allele_loadings$CS1)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot((orderedPC1$CS1)^2, names.arg = row.names(orderedPC1), cex.names = 0.5, las = 2)

#### Getting allele loadings PC2####
orderedPC2<- allele_loadings[order((allele_loadings$CS2)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot((orderedPC2$CS2)^2, names.arg = row.names(orderedPC2), cex.names = 0.5, las = 2)


## PC2 Vs PC3


pca <- qplot(Axis2,Axis3, data = Micro_PCA, color = pop(Microsats_all), cex = 2)
pca + annotate("text", x = Micro_PCA$Axis2, y = Micro_PCA$Axis3, label = row.names(Micro_PCA), cex = 3)

## PC3 is a bit more difficult to interpret, seems to be driven by the gibel goldfish but also PRO samples fall strangely.
## - I think this is probably driven by the GF29 locus and the 191 allele found otherwise only in gibel. 

#### Getting allele loadings PC3 ####
orderedPC3<- allele_loadings[order((allele_loadings$CS3)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot((orderedPC3$CS3)^2, names.arg = row.names(orderedPC3), cex.names = 0.5, las = 2)


## PC3 Vs PC4

### Become less informative after this . . . 
pca <- qplot(Axis3,Axis4, data = Micro_PCA, color = pop(Microsats_all), cex = 2)
pca + annotate("text", x = Micro_PCA$Axis3, y = Micro_PCA$Axis4, label = row.names(Micro_PCA), cex = 3)

## Axis 4 again bringing out some variation between goldfish and gibel

#### Getting allele loadings PC4 ####
orderedPC4<- allele_loadings[order((allele_loadings$CS4)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot((orderedPC4$CS4)^2, names.arg = row.names(orderedPC4), cex.names = 0.5, las = 2)

## PC4 Vs PC5

pca <- qplot(Axis4,Axis5, data = Micro_PCA, color = pop(Microsats_all), cex = 2)
pca + annotate("text", x = Micro_PCA$Axis4, y = Micro_PCA$Axis5, label = row.names(Micro_PCA), cex = 3)

## PC5 gets the within crucian variation, with PRO and STEC being the extremes.

#### Getting allele loadings PC4 ####
orderedPC5<- allele_loadings[order((allele_loadings$CS5)^2, decreasing = T),] ## ordered and squaring to get rid of negative values
barplot((orderedPC5$CS5)^2, names.arg = row.names(orderedPC5), cex.names = 0.5, las = 2)


