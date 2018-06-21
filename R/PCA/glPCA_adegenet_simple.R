#install.packages("adegenet") ## install package, only need to do this once, then you can comment out this line 

library(adegenet) ## load package into this session

setwd("~/Data/RADseq/Rarvalis/Final_populations_outs/") ## wherever you put the input file.


## Reading in the data -------------------------------------------------------------------------------


?read.PLINK ## Get help on the read.PLINK function if you want to take a look
# Basically it makes a "genlight" object which adagenet works with

alldata <- read.PLINK('batch_1.plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=1, saveNbAlleles=T) ## If you are on linus you can change n.cores. But >1 core is not supported on windows.

## A few checks on the data
summary(alldata)
pop(alldata)
indNames(alldata)
nInd(alldata)
ploidy(alldata)

### A useful plot to check the data

plot(alldata, posi = 'topleft', cex = 0.5)

## Also, you can look at the data as a matrix, and it can be subsetted the same way using [indx.row, idx.col] where rows are individuals and columns are loci

#head(as.matrix(alldata)) 


### If you want to add new population names #### -----------------------------------------------------

#pops <- read.delim('pops.txt', header = F) ## Just a single column text file
#samples <- read.delim("sample_names.txt", header = F) ## again ust a single column text file

#pop(alldata) <- pops$V1  ## Assign to the relevant slot in the 
#indNames(alldata) <- samples$V1

##-----------------------------------------------------------------------------------------------------

## PCA analyses ---------------------------------------------------------------------------------------

?glPca ## for help on PCA function. There are some important methodological considerations in here so good to have a look in the adegenet genomics tutorial too. 

pca1 <- glPca(alldata, parallel = FALSE) ## again if on linux, you can set "parallel" to TRUE to use multi cores

scatter(pca1, posi = "topright") ## quick visualisation of first 2 PCs

### Better visualisation

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0,col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2,posi='topright', inset=.05,ratio=.3)



