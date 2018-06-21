
### SETTING UP THE WORKSPACE ## ----------------------------------------------------------

install.packages("polysat") # <- Dunja, this package is already installed on you PC so you don't need to run this line
library(polysat)  # Load this package into the current workspace
install.packages("combinat")  # <- Dunja, you still need to install this package on your PC.
library(combinat)

par(mar = c(4,4,4,4)) # Set display parameters for plotting figures (optional step)

setwd("~/Dropbox/PhD/Dans_PhD_Shared/Triploid_gibel_hybrids/data/") ## Set the working directory. Dunja, change this to the folder that contains your data.

?read.Structure  ## Use this syntax to get help on any function

### DOING THE ANALYSES ### --------------------------------------------------------------

## First, read in from structure file - need to specify ploidy- ie number of lines in the data file per indiv.

data <- read.Structure("1.12.15_ Microsat_Finland_DKL_korriegiert3_NP.txt", ploidy = 4) ## If your file is the same name and you have set the working directory correctly above, the you don't need to change this line.

summary(data) ## Check that the data has been read in correctly optional but advised

## then, after reading in the data, run the estimate ploidy function and manually change the population ploidies below
## Change ploidy of S pop from 4 to 3
## Change ploidy of XO264 from 3 to 2
## Change ploidy of HV1019 from 3 to 2

data <- estimatePloidy(data)

Ploidies(data) ## This is a function just to access the "Ploidies" slot of the polysat object that we have here called "data". Use this to check that the manual changes have worked Ok.

# If we want to use the "Bruvo" distance index in the next step (which is default), we need to add the repeat length for each locus in as below . . . I imagine doing this would be more accurate
Usatnts(data) <- c(2,2,2,2,2,2)

## Now we need to calculate the genetic distances between samples
testmat <- meandistance.matrix(data) ## This should print a load of stuff to the console and may run for a minute or so. 

## And now we can do a principle coordinate analysis on this
pca2 <- cmdscale(testmat, eig = TRUE, k = 5)


## Plotting the PCA ---------------------------------------------------------------------------------------------------


## Note - I have fixed the problem with the colouring, so it will now colour the points by 
## whatever population labels are in the population column of the input file. So to change the colour scheme
## it is easiest to change this line and re-run the script. 

## Here we make an object called "palette" that contains the colours we want to use for plotting. If you want to change these colours you can change this line. But google the colours that you can use in R, there are a lot!!
palette = c("black", "blue", "darkblue", "green", "turquoise", "magenta", "green", "red", "yellow", "purple") ## add as many colours as there are populations

plot(pca2$points[,1], pca2$points[,2], pch = 21, col="black", bg = as.factor(PopInfo(data)), main = "PCA with Bruvo distance", cex = 2, xlab = "PC1", ylab = "PC2") ## Looking at 1st and second PC
text(pca2$points[,1], pca2$points[,2], Samples(data), cex = 0.7) ## Add the sample labels. To change the size, change the value of the cex argument
## My interpretation is that PC1 contains the variation between crucian and gibel
## And PC2 contains variation within gibel

## Can look at other combinations of PCs by changing the numbers in the square brackets. E.g. for PC 2 and 3:
plot(pca2$points[,2], pca2$points[,3], pch = 21, col="black", bg = as.factor(PopInfo(data)), main = "PCA with Bruvo distance", cex = 2)
text(pca2$points[,2], pca2$points[,3], Samples(data), cex = 0.7)



## Examining allele frequencies to find diagnostic alleles. -------------------------------------------------------------------------------------

# We can use the "simpleFreq" function to get allele freqs 
# But remember: "For mixed ploidy systems, the simpleFreq function is available, but will be biased toward underestimating common allele frequencies and overestimating rare allele frequencies, which will cause an underestimation of FST" 
# this isn't too much of a probelm for us if we just want to find diagnostic alleles.

## For the analyses below I changed the population names in the input files to separate crucian and gibel into 2 "Populations".
## This allows us to get allele frequencies for each species separately.
## You could do this for trip and tetraploid gibel too for example.
## Note this was done with the old data file, not the ammended one.

cru_sep_data <- read.Structure("gibel_crucian_separate_pop.txt", ploidy = 4) ## read in the data with new population names

cru_sep_data <- estimatePloidy(cru_sep_data) ## Do the manual ploidy edit again
Ploidies(cru_sep_data) ## Check that it worked

freqs <- simpleFreq(cru_sep_data) # Get the frequency information
freqM <- as.matrix(freqs) # convert this data to a matrix object type (requirement for the barplot function below)

## As there are a lot of alleles, it is a lot to plot on one graph, so I have split the data in half
freqM1 <- freqM[, c(1:(ncol(freqM)/2))] 
freqM2 <- freqM[, c(((ncol(freqM)/2)+1):ncol(freqM))]

## Set the colours to use, here I have set everything to grey except crucian, the order of the colours must correspond to the order of the populations in the input file here!
cols= c("grey", "grey", "grey", "grey", "grey","grey", "grey","grey","red","grey") # Set up the colours to be used for the plot. 

par(mfrow = c(2,1)) ## Set the plot area so it has two panels

## Now plot the two halves of the data
barplot(freqM1[,c(2:ncol(freqM1))], beside = T, col = cols[as.factor(row.names(freqM1))], las = 2, space = c(0,3))
barplot(freqM2[,c(2:ncol(freqM2))], beside = T, col = cols[as.factor(row.names(freqM2))], las = 2, space = c(0,3))




## Calculating Fst's if you want to ---------------------------------------------------------------------------------------------------------

fsts <- calcFst(freqs)


## Note that apparently we should be dealing with allotetraploid genomes separately, by coding the loci for which allele came from which genome . . 
# this should be relatively easy in our data . . . but right now i don't have time . . . 

